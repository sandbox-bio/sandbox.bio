// Defines main logic for parsing and executing commands, using both coreutils
// for the small utilities and Aioli for the bioinformatics tools.

// Known Limitations:
// * Doesn't support subshells, process substitutions, variables, globbing
// * Ctrl + C doesn't stop running process (e.g. `sleep 2 & sleep 5` + ^C does nothing)
// * tail doesn't support `tail -n +3` format (but `head -n-2` and `head/tail -2` supported)


// Imports
import { readable } from "svelte/store";
import parse from "shell-parse";         // Transforms a bash command into an AST
import columnify from "columnify";       // Prettify columnar output
import prettyBytes from "pretty-bytes";  // Prettify number of bytes
import minimist from "minimist";         // Parse CLI arguments
import Aioli from "@biowasm/aioli";

// State
let _aioli = {};   // Aioli object
let _fs = {};      // Aioli filesystem object
let _jobs = 0;     // Number of jobs running in background
let _pid = 10000;  // Current pid


// =============================================================================
// Initialize
// =============================================================================

// Initialize Aioli / prepare bioinformatics tools
// Format: { tools: [], files: [] }
async function init(config={})
{
	_aioli = await new Aioli(config.tools, { env: "stg", debug: false });
	_fs = _aioli.tools[1].module.FS;

	// Pre-load files onto the main folder
	for(let file of config.files || [])
		_fs.writeFile(file.name, file.contents);
}


// =============================================================================
// Execute a command
// =============================================================================

// Run a command
async function exec(cmd, callback)
{
	// console.log("[exec]", cmd, callback);

	// -------------------------------------------------------------------------
	// Input validation
	// -------------------------------------------------------------------------
	if(cmd.type == "subshell")
		throw "Error: subshells are not supported.";
	if(cmd.type == "variableAssignment")
		throw "Error: Bash variables are not supported.";
	if(cmd.body)
		throw `Error: ${cmd.type} is not supported.`;
	if(cmd.args && cmd.args.find(a => a.type == "glob"))
		throw "Error: globbing (e.g. `ls *.txt`) is not supported.";
	if(cmd.args && cmd.args.find(a => a.type == "processSubstitution"))
		throw "Error: process substitution (e.g. `cmd1 <(cmd2)`) is not supported.";

	// If string, convert it to an AST
	if(typeof cmd === "string")
		return await exec(parse(cmd), callback);

	// -------------------------------------------------------------------------
	// Parse commands in AST sequentially
	// -------------------------------------------------------------------------
	if(Array.isArray(cmd))
	{
		let output = "";
		for(let command of cmd)
		{
			// If it's meant to be asynchronous, don't await around; call callback on its own time
			if(command.control == "&")
			{
				const summary = `[${_jobs++}] ${_pid++} `;
				exec(command, callback).then(out => {
					callback(out);
					callback(summary + "done");
					_jobs--;
				});
				callback(summary + "launched");
				continue;
			}

			// Otherwise, run the steps one after the other
			// Note: `&&` and `||` are  handled below
			output += await exec(command, callback);
		}
		return output;
	}

	// -------------------------------------------------------------------------
	// Special "time" command
	// -------------------------------------------------------------------------
	if(cmd.type == "time")
	{
		const timeStart = window.performance.now();
		const output = await exec(cmd.command, callback);
		const timeEnd = window.performance.now();
		callback(`Runtime: ${timeEnd - timeStart}ms`);
		return output;
	}

	// -------------------------------------------------------------------------
	// Process a single command within the AST
	// -------------------------------------------------------------------------
	if(cmd.type == "command")
	{
		// console.log("Command", cmd.command, cmd.args)

		// Interpret the command
		try {
			const tool = cmd.command.value;
			const args = minimist(cmd.args.map(arg => arg.value), minimistConfig[tool]);

			console.log("INTERPRET", tool, args)
			if(tool in coreutils)
				return await coreutils[tool](args);

			// // Otherwise, try running the command with Aioli
			// } else {
			// 	try {
			// 		output = await aioli.exec(cmd);
			// 	} catch (error) {
			// 		output = error;
			// 	}

			// TODO: Handle cmd.next
			// TODO: Handle cmd.redirects

		// If something fails, behave differently depending on e.g. `&&` vs `||`
		} catch (error) {
			// If using `||`, output error but keep going
			if(cmd.control == "||" && cmd.next) {
				callback(error);
				return await exec(cmd.next, callback);
			}
			// Reaches this line if e.g. using `&&`
			throw error;
		}

		return "Unrecognized command:", cmd?.command?.value;
	}

	console.error("Unrecognized command:", cmd);
	throw "Unrecognized command";
}


// -----------------------------------------------------------------------------
// GNU Coreutils JavaScript implementation
// -----------------------------------------------------------------------------

// Utility functions to support coreutils
const utils = {
	// Read entire file contents from virtual file system, and can subset by line numbers
	// e.g. lineRange=`[0]`, `[0,1]`, `[0,5]`
	readFiles: async (paths, lineRange) => {
		if(!Array.isArray(paths))
			paths = [ paths ];

		// For each path, read entire file contents
		const outputs = await Promise.all(paths.map(async path => {
			let output = await _fs.readFile(path, { "encoding": "utf8" });
			// Then subset by lines of interest (not efficient at all, but `FS.read()` fails due to WebWorker issues)
			if(Array.isArray(lineRange))
				output = output.trim().split("\n").slice(...lineRange).join("\n");
			return output;
		}));
		return outputs.join("");
	},
};

// Custom minimist configs for certain coreutils
const minimistConfig = {
	wc: { boolean: ["l", "c", "w"] }
};

// Actual coreutils (only list functions here that should be directly callable by the user from the CLI)
const coreutils = {
	// sleep [1]
	sleep: args => {
		const duration = parseInt(args._[0]) || 1;
		return new Promise((resolve, reject) => {
			setTimeout(() => resolve(""), duration * 1000);
		});
	},

	// cat file1 [file2 ... fileN]
	cat: args => utils.readFiles(args._),

	// head [-n 10|-10]
	head: (args, tail=false) => {
		// Support older `head -5 myfile` format ==> args={ 5: myfile, _: [] }
		const nOld = Object.keys(args).find(arg => parseInt(arg));
		if(nOld)
			args = { n: nOld, _: [ args[nOld] ] };

		// Get first n lines
		const n = args.n || 10;
		return utils.readFiles(args._, tail ? [-n] : [0, n]);
	},
	// tail [-n 10|-10]
	tail: args => coreutils.head(args, true),

	// wc [-l|-w|-c] myfile
	wc: async args => {
		// Count number of lines/words/chars
		const output = await utils.readFiles(args._);
		const nbLines = output.trim().split("\n").length;
		const nbWords = output.trim().split(/\s+/).length;
		const nbChars = output.length;

		// Return result
		if(args.l) return nbLines;
		else if(args.w) return nbWords;
		else if(args.c) return nbChars;
		else return `${nbLines}\t${nbWords}\t${nbChars}`;
	},

	// Small file system utilities
	cd: args => _fs.chdir(args._[0]) && "",
	mv: args => _fs.rename(args._[0], args._[1]) && "",
	rm: args => Promise.all(args._.map(async arg => await _fs.unlink(arg))),
	pwd: args => _fs.cwd(),
	echo: args => args._.join(" "),
	mkdir: args => Promise.all(args._.map(async arg => await _fs.mkdir(arg))),
	rmdir: args => Promise.all(args._.map(async arg => await _fs.rmdir(arg))),
};


// -----------------------------------------------------------------------------
// Export as readable store
// -----------------------------------------------------------------------------

export const CLI = readable({
	init,
	exec
});
