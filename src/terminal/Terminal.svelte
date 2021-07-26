<script>
// Imports
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";
import XTerm from "./XTerm.svelte";


// =============================================================================
// State
// =============================================================================

// Xterm.js
let ready = false;
let term;

// Aioli
let CLI;


// =============================================================================
// On load
// =============================================================================
onMount(async () => {
	// Initialize Aioli
	CLI = await new Aioli(["samtools/1.10", "seqtk/1.2"], { env: "stg", debug: true });

	// Now we can allow user input
	ready = true;
})


// =============================================================================
// Run commands
// =============================================================================

async function exec(cmd, callback)
{
	let out = "";

	// -------------------------------------------------------------------------
	// Support basic coreutils commands
	// -------------------------------------------------------------------------

	// ls
	// FIXME: what if give options to ls?
	if(cmd == "ls" || cmd.startsWith("ls "))
	{
		const folder = cmd.trim().split(" ")[1];
		const output = await CLI.ls(folder || ".");
		console.log(output)
		if(output.mode)
			out = `${output.size}\t${output.mtime}\t${folder}`;
		else if(output)
			out = output.join("\n");
		else
			out = `${folder}: No such file or directory`;

	// Otherwise, try running the command with Aioli
	} else {
		try {
			out = await CLI.exec(cmd.trim());
		} catch (error) {
			out = error;
		}
	}

	term.writeln(out);

	callback();
}
</script>

<XTerm
	ready={ready}
	bind:term={term}
	on:exec={event => exec(event.detail.cmd, event.detail.callback)}
/>
