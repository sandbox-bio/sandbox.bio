<script>
import { onMount, createEventDispatcher } from "svelte";
import { Terminal } from "xterm";
import { WebLinksAddon } from "xterm-addon-web-links";
import { SerializeAddon } from "xterm-addon-serialize";
import { FitAddon } from "xterm-addon-fit";
import LocalEchoController from "local-echo";
import "xterm/css/xterm.css";

// Imports
import { CoreUtils } from "./coreutils";

// Constants
const ANSI_CLEAR = "\x1bc";                // Clear terminal
const dispatch = createEventDispatcher();  // Dispatch for sending "exec" messages to parent component

// Autocomplete subcommands
const AUTOCOMPLETE = {
	samtools: () => ["view", "sort", "index", "idxstats"],
	bedtools2: () => ["intersect", "merge", "complement", "bamtobed"],
	ls: async args => {
		const pathSearch = args[0];                                               // /samtools/examples/toy
		const pathBase = pathSearch.substring(0, pathSearch.lastIndexOf("/")+1);  // /samtools/examples/
		const files = await CoreUtils.ls([pathBase], true);
		return files.map(f => `${pathBase}${f.name}`);
	}
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

	// Disable local-echo's handling of tabs for autocompletion. See handleAutocomplete() for more info.
	addons.echo.handleData_ = addons.echo.handleData;
	addons.echo.handleData = (data) => {
		if(data == "\t")
			return;
		return addons.echo.handleData_(data);
	}

	// Register handlers
	term.onKey(handleShortcuts);
	term.onData(handleAutocomplete);

	// Prepare UI but don't allow input yet
	term.open(divTerminal);
	addons.fit.fit();
});


// =============================================================================
// xterm.js 
// =============================================================================

// Get user input 
function input(toPrint)
{
	if(toPrint)
		term.write(toPrint);
	term.focus();
	addons.echo.read("$ ")
		.then(exec)
		.catch(console.error);
}

// Execute command
async function exec(cmd)
{
	console.log(cmd);

	// -------------------------------------------------------------------------
	// Input validation
	// -------------------------------------------------------------------------
	if(!cmd || !(typeof cmd === "string"))
		return input();
	cmd = cmd.trim();

	// Handle terminal-related commands
	if(cmd == "clear")
		return input(ANSI_CLEAR);

	// -------------------------------------------------------------------------
	// Support basic file redirection (cmd > file) and piping (cmd | cmd2 | cmd3)
	// -------------------------------------------------------------------------

	// Prepare defaults
	let options = {
		cmd: cmd,
		file: null,  // null == stdout, otherwise save to a file path
		pipes: []
	};

	// Redirections: assume ">" is not used in string arguments
	const redirections = cmd.split(">").map(d => d.trim());
	if(redirections.length > 2)
		return input("Unsupported command: Only support one '>' redirection.\n");
	else if(redirections.length == 2) {
		options.cmd = redirections[0];   // e.g. samtools view -q20 toy.sam
		options.file = redirections[1];  // e.g. toy.filtered.sam
	}

	// -------------------------------------------------------------------------
	// Execute command
	// -------------------------------------------------------------------------
	// Send message to parent component asking to execute the command.
	// The callback will output the result and ask for the next input.
	dispatch("exec", {
		cmd: options.cmd,
		callback: async out => {
			let stdout = "";

			// Do we want to save this to a file?
			if(options.file) {
				await CoreUtils.FS.writeFile(options.file, out);
				stdout = "";
			// Or just to stdout
			} else {
				stdout = "" + out;          // Convert to string if it's not already
				if(!stdout.endsWith("\n"))  // Append \n if needed
					stdout += "\n";
			}

			// Output to stdout and ask for next input
			input(stdout);
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
}

// Autocomplete for subcommands. For some reason, it doesn't seem to work with local-echo. When I did "samtools <TAB>",
// it would output "samtools samtools <TAB>" along with the subcommands :(
async function handleAutocomplete(data)
{
	// Parse user input (only care about tabbing to autocomplete)
	if(data != "\t")
		return;
	const input = addons.echo._input;
	const prgm = input.split(" ")[0].trim();  // ls /samtools/examples/toy
	const args = input.split(" ").slice(1);   // ["/samtools/examples/toy"]
	let cacheAutocomplete = [];

	// Turn off local-echo so we can override it's default behavior
	addons.echo.detach();

	// Autocomplete main commands
	// e.g. "<TAB>" --> "samtools   bedtools2"
	if(args.length == 0)
		cacheAutocomplete = Object.keys(AUTOCOMPLETE).filter(d => d.startsWith(prgm));

	// Autocomplete subcommands
	// e.g. "samtools " --> "view   sort   index  idxstats"
	// e.g. "samtools i" --> "index   idxstats"
	else if(args.length < 2 && prgm in AUTOCOMPLETE)
		cacheAutocomplete = (await AUTOCOMPLETE[prgm](args)).filter(d => d.startsWith(args[0]));

	// Process autocomplete
	// If nothing to autocomplete, just add a space
	if(cacheAutocomplete.length == 0) {
		if(!input.endsWith(" "))
			addons.echo.handleCursorInsert(" ");
	}
	// If only one autocomplete option, then autocomplete it!
	else if(cacheAutocomplete.length == 1) {
		const remainingFragment = cacheAutocomplete[0].slice(args.length == 0 ? prgm.length : args[0].length);
		const extraSpacing = (prgm == "ls" && args.length != 0) ? "" : " ";
		addons.echo.handleCursorInsert(remainingFragment + extraSpacing);
	// Otherwise, output all candidates
	} else {
		addons.echo.printAndRestartPrompt(() => addons.echo.printWide(cacheAutocomplete));
	}

	// Re-enable local-echo
	addons.echo.attach();
}
</script>

<div bind:this={divTerminal} style="opacity: { ready ? 1 : 0.6 }"></div>
