<script>
// Imports
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";

import { CoreUtils } from "./coreutils";
import XTerm from "./XTerm.svelte";


// =============================================================================
// State
// =============================================================================

let CLI;            // Aioli object
let ready = false;  // Ready for user input?


// =============================================================================
// On load
// =============================================================================
onMount(async () => {
	// Initialize Aioli
	CLI = await new Aioli(["samtools/1.10", "seqtk/1.2"], { env: "stg", debug: true });
	CoreUtils.CLI = CLI;

	// Now we can allow user input
	ready = true;
});


// =============================================================================
// Run commands
// =============================================================================

async function exec(cmd, callback)
{
	// -------------------------------------------------------------------------
	// Break down command
	// -------------------------------------------------------------------------
	const prgm = cmd.split(" ")[0];
	const args = cmd.split(" ").slice(1);

	// -------------------------------------------------------------------------
	// Support basic coreutils commands
	// -------------------------------------------------------------------------
	if(prgm in CoreUtils) {
		callback(await CoreUtils[prgm](args));
		return;
	}

	// // FIXME: what if give options to ls?
	// let out = "";
	// if(cmd == "ls" || cmd.startsWith("ls "))
	// {
	// 	const folder = cmd.split(" ")[1];
	// 	const output = await CLI.ls(folder || ".");
	// 	console.log(output)
	// 	if(output.mode)
	// 		out = `${output.size}\t${output.mtime}\t${folder}`;
	// 	else if(output)
	// 		out = output.join("\n");
	// 	else
	// 		out = `${folder}: No such file or directory`;

	// Otherwise, try running the command with Aioli
	let out;
	try {
		out = await CLI.exec(cmd);
	} catch (error) {
		out = error;
	}

	callback(out);
}
</script>

<XTerm
	ready={ready}
	on:exec={event => exec(event.detail.cmd, event.detail.callback)}
/>
