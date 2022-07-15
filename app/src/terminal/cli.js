// Defines main logic for parsing and executing commands, using both coreutils
// for the small utilities and Aioli for the bioinformatics tools.
//
// Known Limitations:
// * Doesn't support subshells
// * Ctrl + C doesn't stop running process (e.g. `sleep 2 & sleep 5` + ^C does nothing)

// Imports
import { get, readable } from "svelte/store";
import parse from "shell-parse";         // Transforms a bash command into an AST
import columnify from "columnify";       // Prettify columnar output
import prettyBytes from "pretty-bytes";  // Prettify number of bytes
import minimist from "minimist";         // Parse CLI arguments
import ansiRegex from "ansi-regex";      // Regex for all ANSI codes
import localforage from "localforage";
import Aioli from "@biowasm/aioli";
import { xtermAddons } from "./xterm";
import { env, getLocalForageKey, MAX_FILE_SIZE_TO_CACHE } from "../stores/config";

// State
let _aioli = {};   // Aioli object
let _fs = {};      // Aioli filesystem object
let _jobs = 0;     // Number of jobs running in background
let _pid = 10000;  // Current pid
let _wd = null;    // Track the last folder we were in to support "cd -"

// Define $env for convenience (since not in a .svelte file)
let $env = {};
env.subscribe(d => $env = d);

const DIR_ROOT = "/shared/data";
const DIR_TUTORIALS = `${DIR_ROOT}/tutorials`;
const ESCAPED_EQUAL_SIGN = "___SANDBOXBIO_EQUAL_SIGN___";


// =============================================================================
// Initialize
// =============================================================================

// Initialize Aioli / prepare bioinformatics tools
// Format: { tools: [], files: [], pwd: "bla" }
async function init(config={})
{
	// Initialize
	_aioli = await new Aioli(config.tools, {
		env: ["localhost", "dev.sandbox.bio"].includes(window.location.hostname) ? "stg" : "prd",
		// debug: window.location.hostname == "localhost",
		printInterleaved: false
	});
	_fs = _aioli.tools[1].module.FS;

	// Pre-load files onto the main folder (this happens *before* we load the filesystem state
	// so at this point, /shared/data is empty!)
	if(config.files?.length > 0)
	{
		console.log("Preloading tutorial files...");

		// Setup folders
		const pathDest = `${DIR_TUTORIALS}/${config.pwd || ""}`;
		await exec(`mkdir ${DIR_TUTORIALS} ${pathDest}`);
		await exec(`cd ${pathDest}`);

		// Loop through files (e.g. "data/terminal-basics/orders.tsv")
		for(let file of config.files)
		{
			// Get filename without path (e.g. "orders.tsv")
			const filename = file.split("/").pop();

			// Mount URL (don't need to check whether it exists first since this is empty)
			const url = `${window.location.origin}/${file}`;
			const [path] = await _aioli.mount([ url ]);  // e.g. ["/shared/data/localhost:5000-data-terminal-basics-orders.tsv"]
			// Rename Aioli-mounted URLs that are automatically given long names
			await exec(`mv ${path} ${filename}`);				
		}
	}
}


// =============================================================================
// Custom transformations to apply to CLI commands (for now just aliases)
// =============================================================================

async function transform(cmd)
{
	const tool = cmd.command.value;

	// Handle aliases
	if(tool == "bowtie2")
		cmd.command.value = "bowtie2-align-s";
	else if(tool == "awk")
		cmd.command.value = "gawk";
	else if(tool == "ll") {
		cmd.command.value = "ls";
		cmd.args.push({type: "literal", value: "-l"});
	}

	return cmd;
}


// =============================================================================
// Execute a command
// =============================================================================

// Run a command
async function exec(cmd, callback=console.warn)
{
	// console.log("[exec]", cmd);

	// -------------------------------------------------------------------------
	// Input validation
	// -------------------------------------------------------------------------
	if(cmd.type == "subshell")
		throw "Error: subshells are not supported.";
	if(cmd.body)
		throw `Error: ${cmd.type} is not supported.`;

	// If string, convert it to an AST
	if(typeof cmd === "string") {
		// FIXME: This is a hack to support `awk -v abc=123`, which otherwise breaks
		const awkVarMatch = cmd.matchAll("-v [A-Za-z0-9]+\=[A-Za-z0-9]+");
		if(cmd.includes("awk") && awkVarMatch) {
			for(let match of awkVarMatch)
				cmd = cmd.replace(match[0], match[0].replace("=", ESCAPED_EQUAL_SIGN));
		}

		try {
			return await exec(parse(cmd), callback);
		} catch (error) {
			if(error.name == "SyntaxError")
				throw "Unrecognized command";
			throw error;
		}
	}

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
					callback(summary + "done\n");
					_jobs--;
				});
				callback(summary + "launched\n");
				continue;
			}

			// Otherwise, run the steps one after the other
			// Note: `&&` and `||` are  handled below
			output += await exec(command, callback);
		}
		return output;
	}

	// -------------------------------------------------------------------------
	// Special commands
	// -------------------------------------------------------------------------
	// Estimate command runtime
	if(cmd.type == "time") {
		const timeStart = window.performance.now();
		const output = await exec(cmd.command, callback);
		const timeEnd = window.performance.now();
		callback(`Runtime: ${Math.round((timeEnd - timeStart) * 100)/100}ms\n`);
		return output;
	}

	// Assign a value to a variable
	if(cmd.type == "variableAssignment") {
		const name = cmd.name;
		const value = cmd.value == null ? "" : await utils.getValue(cmd.value);
		env.set({ ...$env, [name]: value });
		return "";
	}

	// -------------------------------------------------------------------------
	// Process a single command within the AST
	// -------------------------------------------------------------------------
	if(cmd.type == "command")
	{
		// console.log("[command]", cmd);

		// Interpret the command
		if(cmd.command == "[[")
			throw "Error: Test expressions not supported";
		try {
			let output;

			// Handle custom transforms on the command
			cmd = await transform(cmd);

			// Parse args
			const tool = cmd.command.value.trim();  // trim removes \n that get introduced if use \
			const argsRaw = (await Promise.all(cmd.args.map(utils.getValue))).flat();
			const args = minimist(argsRaw);

			// If it's a JS-rewritten coreutils, run it
			if(tool in coreutils) {
				output = await coreutils[tool](args);
				if(typeof output === "string" && !output.endsWith("\n"))
					output += "\n";
			// Otherwise, try running the command with Aioli
			} else {
				const argsClean = argsRaw.map(d => d.replaceAll(ESCAPED_EQUAL_SIGN, "="));
				const outputAioli = await _aioli.exec(tool, argsClean);
				const redirect = (cmd.redirects || [])[0];

				// Handle "2>&1" redirection (doesn't yet support "redirectFd" such as 2>somefile)
				if(redirect?.type === "duplicateFd" && redirect?.srcFd === 2 && redirect?.destFd === 1) {
					// Example: "fastp 2>&1 | grep json": at this point we ran fastp and the usage is stored in `outputAioli.stderr`.
					// Instead of outputting it to screen, we add it to the output and using .shift(), remove the 2>&1 operation so
					// we can process the next redirection, if any.
					output = outputAioli.stderr + outputAioli.stdout;
					cmd.redirects.shift();

				// Otherwise, output stderr to screen immediately
				} else {
					callback(outputAioli.stderr);
					// Either output the stdout or pass it along with the pipe
					output = outputAioli.stdout;
				}
			}

			// -----------------------------------------------------------------
			// Handle redirects, e.g. `... | head`, `... > test.txt`
			// -----------------------------------------------------------------
			for(let redirect of (cmd.redirects || []))
			{
				// Handle pipes
				if(redirect.type == "pipe") {
					// Save `output` to a temp file
					const pathTmpFile = await coreutils.mktemp();
					await utils.writeFile(pathTmpFile, output);
					// Run commands with appending arg called temp file, unless
					// user specified we should use stdin via the argument "-".
					const argStdinIndex = redirect.command.args.findIndex(arg => arg.value == "-");
					if(argStdinIndex != -1)
						redirect.command.args[argStdinIndex].value = pathTmpFile;
					else
						redirect.command.args.push({ type: "literal", value: pathTmpFile });
					return exec(redirect.command, callback);

				// Handle redirection to a file
				} else if(redirect.type == "redirectFd") {
					const path = await utils.getValue(redirect.filename);
					// Write to file
					if(redirect.op == ">")
						await utils.writeFile(path, "" + output);  // convert to string if it's not already
					// Append to file (create it first if doesn't exist)
					else if(redirect.op == ">>") {
						try {
							await _fs.stat(path);
						} catch (error) {
							await utils.writeFile(path, "");
						}
						const contents = await _fs.readFile(path, { "encoding": "utf8" });
						await utils.writeFile(path, contents + output);
					} else {
						throw `Error: ${redirect.op} redirections are not supported.`;
					}
					return "";

				// Otherwise, don't know what to do
				} else {
					throw `Error: ${redirect.type} redirections are not supported.`;
				}
			}

			// -----------------------------------------------------------------
			// Handle next command, e.g. `echo 123 && echo 456`, `echo 123 ; echo 456`
			// But if we had `echo 123 || echo 456` and the LHS didn't throw, then we're done
			// -----------------------------------------------------------------
			if(cmd.next && cmd.control != "||")
				output += await exec(cmd.next, callback);  // FIXME: remove breakline and trim

			// If no redirections, just return the result
			return output;

		// If something fails, behave differently depending on e.g. `&&` vs `||`
		} catch (error) {
			// If using `||`, output error but keep going
			if(cmd.control == "||" && cmd.next) {
				callback(error);
				return await exec(cmd.next, callback);
			}
			// Reaches this line if e.g. using `&&` or enter unrecognized program name
			if(typeof error === "string")
				error += "\n";
			throw error;
		}
	}

	console.error("Unrecognized command:", cmd);
	throw "Unrecognized command";
}


// =============================================================================
// GNU Coreutils JavaScript implementation (only list functions here that should
// be directly callable by the user from the CLI).
// =============================================================================

const coreutils = {
	// sleep [1]
	sleep: args => {
		const duration = parseInt(args._[0]) || 1;
		return new Promise((resolve, reject) => {
			setTimeout(() => resolve(""), duration * 1000);
		});
	},

	// env
	env: args => Object.keys($env).map(v => `${v}=${$env[v]}`).join("\n"),
	hostname: args => "sandbox",
	uname: args => "sandbox.bio",
	whoami: args => $env?.USER || "guest",
	unset: args => {
		args._.map(v => delete $env[v]);
		env.set($env);
		return "";
	},

	// -------------------------------------------------------------------------
	// File system management
	// -------------------------------------------------------------------------
	mv: args => _fs.rename(args._[0], args._[1]) && "",
	rm: args => Promise.all(args._.map(async arg => await _fs.unlink(arg))) && "",
	pwd: args => _fs.cwd(),
	touch: async args => {
		return Promise.all(args._.map(async path => {
			try {
				// Will throw if file doesn't exist
				await _fs.stat(path);
				await _fs.utime(path, new Date().getTime(), new Date().getTime());
			} catch (error) {
				await utils.writeFile(path, "");
			}	
		})) && "";
	},
	cd: async args => {
		let dir = args._[0];
		// Support cd ~ and cd -
		if(dir == "~")
			dir = $env.HOME;
		else if(dir == "-" && _wd)
			dir = _wd;
		// Change directory
		const dirOld = await coreutils.pwd();
		try {
			await _aioli.cd(dir);
			_wd = dirOld;
		} catch (error) {
			return `${dir}: No such file or directory`;
		}
		return "";
	},
	mkdir: args => Promise.all(args._.map(async arg => {
		try {
			await _aioli.mkdir(arg);
		} catch (error) {
			return `${arg}: Cannot create folder`;
		}
	})) && "",
	rmdir: args => Promise.all(args._.map(async arg => await _fs.rmdir(arg))) && "",
	mktemp: args => {
		const path = `/shared/tmp/tmp${parseInt(Math.random() * 1_000_000)}`;
		_fs.writeFile(path, "");
		return path;
	},
	cp: async args => {
		// Copy file contents if it exists
		let data;
		try {
			// Read source file
			await utils.ls({_: [args._[0]]});
			data = await _fs.readFile(args._[0], { encoding: "binary" });

			// Copy data over
			if(args._[1] == ".")
				args._[1] = args._[0].split("/").pop();
			await utils.writeFile(args._[1], data, { encoding: "binary" });
			return "";
		} catch (error) {
			return error;
		}
	},

	// open <file>: command that loads a file in a new tab (used by tutorials to open html files)
	open: async args => {
		const file = args._[0];
		const type = file.endsWith(".html") ? "text/html" : "text/plain";

		// Create Blob from contents (can't use _aioli.download b/c need to specify type)
		const contents = await _aioli.cat(file);
		const blob = new Blob([ contents ], { type });
		const url = URL.createObjectURL(blob);
		window.open(url);

		return "";
	},

	// download <file>: command that downloads a file
	download: async args => {
		const file = args._[0];

		// Read file as binary, even if it's a plain text file since we're just downloading it
		const contents = await _fs.readFile(file, { encoding: "binary" });

		// Create Blob from contents (can't use _aioli.download b/c need to specify encoding/type).
		// Set type to "application/octet-stream" so that location.href downloads the file.
		const blob = new Blob([ contents ], { type: "application/octet-stream" });
		const url = URL.createObjectURL(blob);

		// Create link element to customize the filename, otherwise it's a UUID.
		const fileLink = document.createElement("a");
		fileLink.href = url;
		fileLink.download = file;
		fileLink.click();

		return "";
	},

	// -------------------------------------------------------------------------
	// curl <http...>
	// -------------------------------------------------------------------------
	curl: async args => {
		const url = args._[0];
		const contents = await fetch(url).then(d => d.text());

		// Output to file
		if(args.o || args.O) {
			if(args.O)
				args.o = url.split("/").pop();
			await utils.writeFile(args.o, contents);
			return "";
		}

		return contents;
	},
};


// =============================================================================
// Utility functions
// =============================================================================

const utils = {
	// Get the value of an argument (recursive if need-be)
	getValue: async arg => {
		// Literal; support ~
		if(arg.type == "literal") {
			if(arg.value === "~" || arg.value?.startsWith("~/"))
				return arg.value.replaceAll("~", $env.HOME);
			return arg.value;
		// Variable
		} else if(arg.type == "variable")
			return $env[arg.name] || "";
		// Variable concatenation, e.g. echo "something $abc $def"
		else if(arg.type == "concatenation")
			return (await Promise.all(arg.pieces.map(utils.getValue))).join("");
		// Command Substitution, e.g. someprgm $(grep "bla" test.txt | wc -l)
		else if(arg.type == "commandSubstitution")
			return await exec(arg.commands);
		// Process substitution, e.g. bedtools -a <(grep "Enhancer" data.bed)
		else if(arg.type == "processSubstitution" && arg.readWrite == "<") {
			const output = await exec(arg.commands);
			const pathTmpFile = await coreutils.mktemp();
			await utils.writeFile(pathTmpFile, output);
			return pathTmpFile;
		// Globbing, e.g. ls c*.b?d
		} else if(arg.type == "glob") {
			// Get files at the base path
			if(arg.value.endsWith("/"))
				arg.value = arg.value.slice(0, -1);
			const pathBase = arg.value.substring(0, arg.value.lastIndexOf("/") + 1) || "";
			let pathPattern = arg.value.replace(pathBase, "");
			if(pathPattern == "*")
				pathPattern = "";
			const files = await utils.ls([ pathBase || "." ], true);

			// Convert bash regex to JS regex; "*" in bash == ".*" in js; "?" in bash == "." in js
			const pattern = pathPattern
				.replaceAll("*", "###__ASTERISK__###")
				.replaceAll("?", "###__QUESTION__###")
				.replace(/[-\/\\^$*+?.()|[\]{}]/g, "\\$&")  // escape non-regex chars (e.g. the "." in the file extension)
				.replaceAll("###__ASTERISK__###", ".*")
				.replaceAll("###__QUESTION__###", ".");
			// If user specifies ls *txt, match both hello.txt and my_txt/
			const re = new RegExp("^" + pattern + "$|" + "^" + pattern + "/$");

			// If find no matches, return original glob value
			const filesMatching = files.filter(f => f.name.match(re)).map(f => `${pathBase}${f.name}`)
			if(filesMatching.length > 0)
				return filesMatching;
			if(pathPattern == "")
				return [pathBase || "."];
			return arg.value;
		}
		else
			throw `Error: ${arg.type} arguments not supported.`;
	},

	// Write file
	writeFile: async (path, contents, opts={ encoding: "utf8" }) => {
		// Delete destination file if it exists (otherwise can get errors if destination = lazyloaded URL)
		try {
			await utils.ls({_: [ path ]});
			await coreutils.rm({_: [ path ]});
		} catch (error) {}  // don't error if the destination doesn't exist

		// Remove ANSI characters from file contents before saving (only if string; could also be Uint8Array for binary files)
		if(typeof contents === "string")
			contents = contents.replace(ansiRegex(), "");

		await _fs.writeFile(path, contents, opts);
	},

	// Mount file(s)
	mount: async (files) => {
		return await _aioli.mount(files);
	},

	// -------------------------------------------------------------------------
	// ls [-l] <file1> <file2>
	// We need this to e.g. store the FS state in localStorage, and for other coreutils above
	// -------------------------------------------------------------------------
	ls: async (args, raw=false) => {
		// Input validation and dedupping
		let paths = args._ || args;
		if(paths.length == 0)
			paths = ["."];
		paths = [...new Set(paths)];

		// Gather stats about files and folders
		let stats = [];
		for(let path of paths)
		{
			try {
				let stat = await _fs.stat(path);
				// If the path is a file, we already have the info we need
				if(!await _fs.isDir(stat.mode))
					stats.push({
						name: path.split("/").pop(),
						size: stat.size,
						date: stat.mtime
					});

				// But if the path is a folder, get stat for each node in the folder
				else
				{
					const files = await _fs.readdir(path);
					for(let f of files)
					{
						if(f == "." || f == "..")
							continue;
						stat = await _fs.stat(`${path}/${f}`);
						const isDir = await _fs.isDir(stat.mode);
						let name = isDir ? f + "/" : f;
						name = (!raw && isDir) ? `\x1b[1;34m${name}\x1b[0m` : name;
						stats.push({
							name: name,
							size: stat.size,
							date: stat.mtime
						});
					}
				}
			} catch (error) {
				// console.error(error);
				throw `${path}: No such file or directory`;
			}
		}

		// If don't want to ouput to screen
		if(raw)
			return stats;

		// Columnify output
		return columnify(stats.map(f => ({
			...f,
			size: prettyBytes(f.size),
			date: f.date.toLocaleString()
		})), {
			showHeaders: false,
			columnSplitter: "\t",
			columns: ["size", "date", "name"]
		});
	}
};


// =============================================================================
// File system / CLI history caching functions
// =============================================================================

const fsSave = async function() {
	// console.log("Saving filesystem state...")
	const filesToCache = await fsTraverse(`${DIR_ROOT}/`);

	// Cache user-created files in a localforage key
	const files = {}, folders = {};
	for(let path of filesToCache) {
		// Don't cache large files (e.g. accidentally create large temp file or mount large local file)
		// If .stat fails, then skip this file
		try {
			const size = (await _fs.stat(path)).size;
			if(size > MAX_FILE_SIZE_TO_CACHE)
				continue;
		} catch (error) {
			console.warn(error);
			continue;			
		}

		// For folders, just need to know they're there
		if(path.endsWith("/"))
			folders[path] = true;
		else
			files[path] = await _fs.readFile(path, { "encoding": "binary" });
	}
	await localforage.setItem(`${getLocalForageKey("fs")}files`, files);
	await localforage.setItem(`${getLocalForageKey("fs")}folders`, folders);
	
	// While we're here, also save command history
	const history = get(xtermAddons)?.echo?.history?.entries || [];
	await localforage.setItem(getLocalForageKey("history"), history);
}

const fsLoad = async function() {
	console.log("Loading filesystem state...")
	const files = await localforage.getItem(`${getLocalForageKey("fs")}files`);
	const folders = await localforage.getItem(`${getLocalForageKey("fs")}folders`);

	// First, create the folders, then the files they contain
	for(let path in folders) {
		try {
			await _fs.stat(path);
		} catch (error) {
			await exec(`mkdir ${path}`);
		}
	}
	for(let path in files)
		await utils.writeFile(path, files[path], { encoding: "binary" });

	// While we're here, also load saved command history
	const history = await localforage.getItem(getLocalForageKey("history")) || [];
	const historyController = get(xtermAddons)?.echo?.history;
	if(historyController) {
		historyController.entries = history;
		historyController.cursor = history.length;
	}
}

// Recursively traverse folder path and get list of all files
async function fsTraverse(path) {
	let paths = [];
	const files = await utils.ls([path], true);
	for(let file of files) {
		const filePath = path + file.name;
		paths.push(filePath);
		if(file.name.endsWith("/"))
			paths = paths.concat(await fsTraverse(filePath));
	}
	return paths;
}


// =============================================================================
// Export CLI as a readable store
// =============================================================================

export const CLI = readable({
	init,
	exec,
	coreutils,
	utils,
	fsSave,
	fsLoad
});
