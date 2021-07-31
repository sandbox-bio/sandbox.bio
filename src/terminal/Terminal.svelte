<script>
// Imports
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";

import { CoreUtils } from "./coreutils";
import XTerm from "./XTerm.svelte";

export let tools = ["samtools/1.10", "bedtools/2.29.2"];
export let files = [];
export let intro;


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
	try {
		CLI = await new Aioli(tools, { env: "stg", debug: false });
		CoreUtils.CLI = CLI;

		// Now we can allow user input
		ready = true;

		// Pre-load files onto the main folder
		for(let file of files)
			await CoreUtils.FS.writeFile(file.name, file.contents);

	} catch (error) {
		console.error("Could not load biowasm modules:", error);
		alert("Could not load the terminal. Please try refreshing the page.");
	}
});


// =============================================================================
// Run commands
// =============================================================================

async function exec(cmd, callback)
{
	let output = "";
	const prgm = cmd.split(" ")[0];
	const args = cmd.split(" ").slice(1);

	// Is this a coreutils command?
	if(prgm in CoreUtils) {
		try {
			output = await CoreUtils[prgm](args);			
		} catch (error) {
			output = error;
		}

	// Otherwise, try running the command with Aioli
	} else {
		try {
			output = await CLI.exec(cmd);
		} catch (error) {
			output = error;
		}
	}

	callback(output);
}
</script>

<XTerm {ready} {intro} on:exec={event => exec(event.detail.cmd, event.detail.callback)} />
