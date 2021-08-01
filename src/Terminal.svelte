<script>
// Imports
import { onMount } from "svelte";
import XTerm from "./XTerm.svelte";
import { CLI } from "terminal/cli";

// Props
export let tools = ["samtools/1.10", "bedtools/2.29.2"];
export let files = [];
export let intro = false;

// State
let ready = false;  // Ready for user input?


// =============================================================================
// On load
// =============================================================================
onMount(async () => {
	// Initialize Aioli
	try {
		await $CLI.init({ tools, files });
		ready = true;
	} catch (error) {
		console.error("Could not load terminal:", error);
		alert("Could not load the terminal. Please try refreshing the page.");
	}
});


// =============================================================================
// Run commands
// =============================================================================

async function exec(cmd, callback) {
	const output = await $CLI.exec(cmd);
	callback(output);
}
</script>

<XTerm {ready} {intro} on:exec={event => exec(event.detail.cmd, event.detail.callback)} />
