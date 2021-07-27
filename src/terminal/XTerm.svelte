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
// Autocomplete subcommands
const AUTOCOMPLETE = {
	samtools: ["view", "sort", "index", "idxstats"],
	bedtools2: ["intersect", "merge", "complement", "bamtobed"]
};


// =============================================================================
// State
// =============================================================================

export let ready = false;  // Whether CLI is ready for user input
let term;                  // XTerm.js object
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
		echo: new LocalEchoController(null, {                  // Echo controller
			historySize: 1000
		}),
		serialize: new SerializeAddon(),                       // Can be used to save state
		links: new WebLinksAddon(),                            // Supports links in the terminal
		fit: new FitAddon()                                    // Makes terminal fit HTML element
	};

	// Attach addons
	for(let addonName in addons)
		term.loadAddon(addons[addonName]);

	// Register handlers
	term.onKey(handleShortcuts);
	term.onData(handleData);
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
	console.log(cmd);

	// Input validation
	if(!cmd || !(typeof cmd === "string")) {
		input();
		return;
	}
	cmd = cmd.trim();

	// Handle terminal-related commands
	if(cmd == "clear") {
		term.write(ANSI_CLEAR);
		input();
		return;
	}

	// Send message to parent component asking to execute the command.
	// The callback will output the result and ask for the next input.
	dispatch("exec", {
		cmd: cmd,
		callback: out => {
			term.writeln("" + out);  // convert to string if it's not
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
	// Ctrl + L = Clear terminal (also clear input in case user had written something on the line)
	if(key.domEvent.ctrlKey && key.domEvent.key == "l") {
		term.write(ANSI_CLEAR);
		addons.echo.setInput("");
	}

	// Ctrl + A = Beginning of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "a")
		addons.echo.setCursor(0);

	// Ctrl + E = End of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "e")
		addons.echo.setCursor(Infinity);

	// Tab = try to autocomplete bioinformatics subcommands. This is part 1.
	// See handleAutocomplete() below for part 2.
	if(key.domEvent.key == "Tab")
	{

// Autocomplete for subcommands. For some reason, it doesn't seem to work with local-echo. When I did "samtools <TAB>",
// it would output "samtools samtools <TAB>" along with the subcommands :(
function handleData(data)
{
	// Parse user input (only care about tabbing to autocomplete)
	if(data != "\t")
		return;
	const input = addons.echo._input;
	const prgm = input.split(" ")[0].trim();
	const args = input.split(" ").slice(1);
	let cacheAutocomplete = [];

	// Autocomplete subcommands: e.g. "samtools <TAB>" will show "view   sort   index  idxstats"
	if(args.length < 2 && prgm in AUTOCOMPLETE)
	{
		// Turn off local-echo so we can override it's default behavior
		addons.echo.detach();

		// e.g. "samtools " --> "view   sort   index  idxstats"
		// e.g. "samtools i" --> "index   idxstats"
		cacheAutocomplete = AUTOCOMPLETE[prgm].filter(d => d.startsWith(args[0]));

		// If only one autocomplete option, then autocomplete it!
		if(cacheAutocomplete.length == 1)
			addons.echo.handleCursorInsert(cacheAutocomplete[0].slice(args[0].length) + " ");
		// Otherwise, output candidates
		else
			addons.echo.printAndRestartPrompt(() => addons.echo.printWide(cacheAutocomplete));

		// Re-enable local-echo
		addons.echo.attach();
	}
}

// Auto-completes common commands
function handleAutocomplete(index, tokens)
{
	// Root autocomplete: show supported tools
	if(index == 0)
		return Object.keys(AUTOCOMPLETE);

	// This generates an error on purpose.
	return;  // this generates an error
}
</script>

<div bind:this={divTerminal} style="opacity: { ready ? 1 : 0.6 }"></div>
