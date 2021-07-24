<script>
import { onMount } from "svelte";
import { Terminal } from "xterm";
import { WebLinksAddon } from "xterm-addon-web-links";
import LocalEchoController from "local-echo";
import "xterm/css/xterm.css";

// =============================================================================
// State
// =============================================================================

// Terminal state
let term;
let termEcho;

// UI
let divTerminal;

// ANSI Constants: <https://en.wikipedia.org/wiki/ANSI_escape_code#CSI_(Control_Sequence_Introducer)_sequences>
const ANSI_ESC = "\x1b";
const ANSI_CLEAR = ANSI_ESC + "c";


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
	term.loadAddon(new WebLinksAddon());

	// Register handlers
	term.onKey(handleShortcuts);
	termEcho.addAutocompleteHandler(handleAutocomplete);

	// Ready for user input
	term.writeln("# sandbox.bio\n");
	term.focus();
	input();
})


// =============================================================================
// Handlers
// =============================================================================

// Keyboard shortcuts
function handleShortcuts(key)
{
	console.log(key);

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


// =============================================================================
// Manage user input
// =============================================================================

async function exec(cmd)
{
	console.log("Command:", cmd);

	// Basic commands
	if(cmd == "clear") {
		term.write(ANSI_CLEAR);
		return;
	}
}

// Get user input
function input()
{
	termEcho.read("$ ")
		.then(exec)
		.catch(console.error)
		.finally(input);
}
</script>

<div bind:this={divTerminal}></div>
