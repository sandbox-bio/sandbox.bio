<script>
// Imports
import { onMount } from "svelte";
import { Terminal } from "xterm";
import { WebLinksAddon } from "xterm-addon-web-links";
import { SerializeAddon } from "xterm-addon-serialize";
import { FitAddon } from "xterm-addon-fit";
import LocalEchoController from "local-echo";
import Aioli from "@biowasm/aioli";
import "xterm/css/xterm.css";


// =============================================================================
// State
// =============================================================================

// Terminal state
let term;
let termEcho;
let termLinks;
let termSerialize;
let termFit;

// UI
let divTerminal;

// ANSI Constants: <https://en.wikipedia.org/wiki/ANSI_escape_code#CSI_(Control_Sequence_Introducer)_sequences>
const ANSI_ESC = "\x1b";
const ANSI_CLEAR = ANSI_ESC + "c";

// CLI
let ready = false;  // If true, ready for user input
let CLI;            // Aioli object


// =============================================================================
// On load
// =============================================================================
onMount(async () => {
	// Initialize terminal
	term = new Terminal({
		convertEol: true,
		cursorBlink: true
	});
	term.open(divTerminal);

	// Attach addons
	term.loadAddon(termEcho = new LocalEchoController());
	term.loadAddon(termSerialize = new SerializeAddon());  // Can be used to save state
	term.loadAddon(termLinks = new WebLinksAddon());       // Supports links in the terminal
	term.loadAddon(termFit = new FitAddon());              // Makes terminal fit HTML element

	// Register handlers
	term.onKey(handleShortcuts);
	termEcho.addAutocompleteHandler(handleAutocomplete);

	// Initialize Aioli
	CLI = await new Aioli(["samtools/1.10", "seqtk/1.2"], { env: "stg", debug: true });
	// Make Terminal ready for user input
	termFit.fit();
	input();
	term.focus();
	ready = true;
})


// =============================================================================
// Manage user input
// =============================================================================

async function exec(cmd)
{
	if(!cmd || !cmd.trim())
		return;

	console.log("Command:", cmd);

	let out = "";

	// -------------------------------------------------------------------------
	// Support basic coreutils commands
	// -------------------------------------------------------------------------

	// Basic commands
	if(cmd == "clear")
		out = ANSI_CLEAR;
	// ls
	// FIXME: what if give options to ls?
	else if(cmd == "ls" || cmd.startsWith("ls "))
	{
		const folder = cmd.trim().split(" ")[1];
		const output = await CLI.ls(folder || ".");
		console.log(output)
		if(output.mode)
			out = `${output.size}\t${output.mtime}\t${folder}`;
		else if(output)
			out = output.join("\n");
		else
			out = `${folder}: No such file or directory`

	// Otherwise, try running the command with Aioli
	} else {
		try {
			out = await CLI.exec(cmd.trim());
		} catch (error) {
			out = error;
		}
	}

	term.writeln(out + "\n");
}

// Get user input
function input()
{
	termEcho.read("$ ")
		.then(exec)
		.catch(console.error)
		.finally(input);
}


// =============================================================================
// xterm.js handlers
// =============================================================================

// Keyboard shortcuts
function handleShortcuts(key)
{
	// Ctrl + L = Clear terminal
	if(key.domEvent.ctrlKey && key.domEvent.key == "l")
		term.write(`${ANSI_CLEAR}$ `);

	// Ctrl + A = Beginning of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "a")
		termEcho.setCursor(0);

	// Ctrl + E = End of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "e")
		termEcho.setCursor(Infinity);
}

// Auto-completes common commands
function handleAutocomplete(index, tokens)
{
	const command = tokens[0];
	console.log(index, tokens);

	// Root autocomplete
	if(index == 0)
		return ["samtools", "bedtools2"];

	// Samtools autocomplete. Need `&& tokens[index]`, otherwise results in "samtools samtools"
	if(index == 1 && command == "samtools" && tokens[index])
		return ["view", "index", "sort"];

	return [];
}
</script>

<div bind:this={divTerminal} style="opacity: { ready ? 1 : 0.6 }"></div>
