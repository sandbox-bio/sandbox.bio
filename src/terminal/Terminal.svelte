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
	CLI = await new Aioli(["samtools/1.10", "bedtools2/2.29.2"], { env: "stg", debug: true });
	CoreUtils.CLI = CLI;

	// Now we can allow user input
	ready = true;
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

<XTerm
	ready={ready}
	on:exec={event => exec(event.detail.cmd, event.detail.callback)}
/>
