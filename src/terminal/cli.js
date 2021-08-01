// =============================================================================
// Defines main logic for parsing and executing commands, using both coreutils
// for the small utilities and Aioli for the bioinformatics tools.
// =============================================================================

// Imports
import { get, readable } from "svelte/store";
import columnify from "columnify";       // Prettify columnar output
import prettyBytes from "pretty-bytes";  // Prettify number of bytes
import minimist from "minimist";         // Parse CLI arguments
import Aioli from "@biowasm/aioli";

import parse from "shell-parse";

// import { coreutils } from "./coreutils";
import { xterm, xtermAddons } from "./xterm";

// State
let aioli = {};  // Aioli object
let FS = {};     // Aioli filesystem object

// Convenience variables since can't use $store outside .svelte file
const $xtermAddons = get(xtermAddons);


// -----------------------------------------------------------------------------
// Initialize
// -----------------------------------------------------------------------------

// Initialize Aioli / prepare bioinformatics tools
// Format: { tools: [], files: [] }
async function init(config={})
{
	aioli = await new Aioli(config.tools, { env: "stg", debug: false });
	FS = aioli.tools[1].module.FS;

	// Pre-load files onto the main folder
	for(let file of config.files || [])
		FS.writeFile(file.name, file.contents);
}


// -----------------------------------------------------------------------------
// Execute a command
// -----------------------------------------------------------------------------

// Run a command
async function exec(cmd)
{
	console.log("[exec]", cmd);

	// If the input is a string, means it's a string command we need to turn into an AST
	if(typeof cmd === "string")
		return await exec(parse(cmd));

	// If the input is an array, assume it's an AST containing a list of commands to run sequentially
	if(Array.isArray(cmd)) {
		let output = "";
		for(let command of cmd)
		{
			// If needs to be asynchronous
			output += await exec(command);
		}
		return output;
	}

	// Otherwise, we're looking at a single command within an AST
	if(cmd.control == "&")
		throw "Error: Asynchronous commands (e.g. `echo 123 &`) are not supported.\n";
	if(cmd.type == "command")
	{
		console.log("Command", cmd.command, cmd.args)
		// cmd.next
		// cmd.redirects

		return "/done";
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
			await aioli.exec("bedtools")
		)
	}
};


// -----------------------------------------------------------------------------
// Small utilities
// -----------------------------------------------------------------------------

//
function setInput(command) {
	$xtermAddons.echo.setInput(command);
	$xtermAddons.echo.handleData("\r");
}


// -----------------------------------------------------------------------------
// Export as readable store
// -----------------------------------------------------------------------------

export const CLI = readable({
	init,
	exec,
	setInput
});
