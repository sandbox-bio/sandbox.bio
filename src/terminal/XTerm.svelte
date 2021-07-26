<script>
import { onMount, createEventDispatcher } from "svelte";
import { Terminal } from "xterm";
import { WebLinksAddon } from "xterm-addon-web-links";
import { SerializeAddon } from "xterm-addon-serialize";
import { FitAddon } from "xterm-addon-fit";
import LocalEchoController from "local-echo";
import "xterm/css/xterm.css";

// Constants
const ANSI_CLEAR = "\x1bc";                // Clear terminal
const dispatch = createEventDispatcher();  // Dispatch for sending "exec" messages to parent component


// =============================================================================
// State
// =============================================================================

export let ready = false;  // Whether CLI is ready for user input
export let term;           // XTerm.js object

let addons;                // XTerm.js add-ons
let divTerminal;           // HTML element where terminal will be drawn

$: if(ready) input();      // Ask for user input once ready


// =============================================================================
// Initialization
// =============================================================================

onMount(() => {
	// Xterm.js
	term = new Terminal({
		convertEol: true,
		cursorBlink: true
	});

	// Xterm.js add-ons
	addons = {
		echo: new LocalEchoController(),  // Echo controller
		serialize: new SerializeAddon(),  // Can be used to save state
		links: new WebLinksAddon(),       // Supports links in the terminal
		fit: new FitAddon()               // Makes terminal fit HTML element
	};

	// Attach addons
	for(let addonName in addons)
		term.loadAddon(addons[addonName]);

	// Register handlers
	term.onKey(handleShortcuts);
	addons.echo.addAutocompleteHandler(handleAutocomplete);

	// Prepare UI but don't allow input yet
	term.open(divTerminal);
	addons.fit.fit();
});


// =============================================================================
// xterm.js 
// =============================================================================

// Get user input 
function input()
{
	term.focus();
	addons.echo.read("$ ")
		.then(exec)
		.catch(console.error);
}

// Execute command
function exec(cmd)
{
	// Input validation
	if(!cmd || !(typeof cmd === "string"))
		return;
	cmd = cmd.trim();

	// Handle terminal-related commands
	if(cmd == "clear") {
		term.write(ANSI_CLEAR);
		return;
	}

	// Send message to parent component asking to execute the command.
	// The callback will output the result and ask for the next input.
	dispatch("exec", {
		cmd: cmd,
		callback: out => {
			term.writeln(out);
			input();
		}
	});
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
		addons.echo.setCursor(0);

	// Ctrl + E = End of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "e")
		addons.echo.setCursor(Infinity);
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
