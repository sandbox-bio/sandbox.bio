// =============================================================================
// Defines main logic for parsing and executing commands, using both coreutils
// for the small utilities and Aioli for the bioinformatics tools.
// =============================================================================

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


// -----------------------------------------------------------------------------
// Initialize
// -----------------------------------------------------------------------------

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


// -----------------------------------------------------------------------------
// Execute a command
// -----------------------------------------------------------------------------

// Run a command
async function exec(cmd, callback)
{
	console.log("[exec]", cmd, callback);

	// Input validation
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

	// If the input is a string, it's a string command we need to turn into an AST
	if(typeof cmd === "string")
		return await exec(parse(cmd), callback);

	// If the input is an array, it's an AST containing a list of commands to run sequentially
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

	// Otherwise, we're looking at a single command within an AST
	if(cmd.type == "command")
	{
		console.log("Command", cmd.command, cmd.args)
		// cmd.next
		// cmd.redirects

		// Support `&&` and `||`
		try {
			d
		} catch (error) {
			// If using `||`, output error but keep going
			if(cmd.control == "||" && cmd.next) {
				callback(error);
				return await exec(cmd.next, callback);
			}
			// Reaches this line if e.g. using `&&`
			throw error;
		}

		return "/done";
	}

	// If user wants to track runtime
	if(cmd.type == "time")
	{
		const timeStart = window.performance.now();
		const output = await exec(cmd.command, callback);
		const timeEnd = window.performance.now();
		callback(`Runtime: ${timeEnd - timeStart}ms`);
		return output;
	}

	console.error("Unrecognized command:", cmd);
	throw "Unrecognized command";

	// console.log(minimist(cmd.split(" "), {
	// 	unknown: d => {
	// 		console.log("unknown", d)
	// 	}
	// }))
	// console.log(await coreutils.head("yep"))

	// let output = "";
	// const prgm = cmd.split(" ")[0];
	// const args = cmd.split(" ").slice(1);

	// // Is this a coreutils command?
	// if(prgm in CoreUtils) {
	// 	try {
	// 		output = await CoreUtils[prgm](args);			
	// 	} catch (error) {
	// 		output = error;
	// 	}

	// // Otherwise, try running the command with Aioli
	// } else {
	// 	try {
	// 		output = await aioli.exec(cmd);
	// 	} catch (error) {
	// 		output = error;
	// 	}
	// }

	// return output;
}



// -----------------------------------------------------------------------------
// GNU Coreutils JavaScript implementation
// -----------------------------------------------------------------------------

const coreutils = {
	head: async (args) => {
		console.log(
			await _aioli.exec("bedtools")
		)
	}
};


// -----------------------------------------------------------------------------
// Export as readable store
// -----------------------------------------------------------------------------

export const CLI = readable({
	init,
	exec
});
