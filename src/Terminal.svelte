<script>
import { onMount, createEventDispatcher } from "svelte";
import { watchResize } from "svelte-watch-resize";
import { Spinner } from "sveltestrap";
import "xterm/css/xterm.css";

// Imports
import { xterm, xtermAddons } from "terminal/xterm";
import { CLI } from "terminal/cli";

// Constants
const ANSI_CLEAR = "\x1bc";
const TOOLS_DEFAULT = ["samtools/1.10", "bedtools/2.29.2", "bowtie2/bowtie2-align-s/2.4.2"];

// Autocomplete subcommands
const AUTOCOMPLETE = {
	samtools: () => ["view", "sort", "depth", "index", "idxstats", "flags", "flagstats"],
	bedtools: () => ["intersect", "merge", "complement", "genomecov", "jaccard", "makewindows", "flank"],
	bowtie2: () => [],
	ls: async args => {
		let pathSearch = args[0];                                                      // /samtools/examples/toy
		let pathBase = pathSearch.substring(0, pathSearch.lastIndexOf("/")+1) || ".";  // /samtools/examples/
		const files = await $CLI.coreutils.ls([pathBase], true);
		if(pathBase == ".")
			pathBase = "";
		return files.map(f => `${pathBase}${f.name}`);
	},
	cat: () => [],
	head: () => [],
	tail: () => [],
	wc: () => [],
	pwd: () => [],
	cd: () => [],
	echo: () => [],
	mv: () => [],
	rm: () => [],
	mkdir: () => [],
	rmdir: () => [],
	env: () => []
};


// =============================================================================
// State
// =============================================================================

export let ready = false;                  // Whether CLI is ready for user input
export let intro = "";                     // Intro string to display on Terminal once ready (optional)
export let init = "";                      // Command to run to initialize the environment (optional)
export let files = [];                     // Files to preload on the filesystem
export let tools = TOOLS_DEFAULT;          // Aioli tools to load

let divTerminal;                           // HTML element where terminal will be drawn
$: if(ready) input();                      // Ask for user input once ready

const dispatch = createEventDispatcher();  // Send info to parent component when cmd is done


// =============================================================================
// Initialization
// =============================================================================

onMount(async () => {
	// Register handlers
	$xterm.onKey(handleShortcuts);
	$xterm.onData(handleAutocomplete);

	if(intro)
		$xterm.writeln(intro);

	// Prepare UI but don't allow input yet
	$xterm.open(divTerminal);

	// Initialize Aioli
	try {
		await $CLI.init({ tools, files });
		if(init)
			await $CLI.exec(init);
		ready = true;
	} catch (error) {
		console.error("Could not load terminal:", error);
		alert("Could not load the terminal. Please try refreshing the page.");
	}
});

// Resize xterm when the window size changes
function handleResize() {
	$xtermAddons.fit.fit();
}

// =============================================================================
// xterm.js 
// =============================================================================

// Get user input 
function input(toPrint)
{
	if(toPrint)
		$xterm.writeln(toPrint);
	$xterm.focus();
	$xtermAddons.echo.read("$ ")
		.then(exec)
		.catch(console.error);
}

// Execute command
async function exec(cmd)
{
	// Handle special `clear` command
	if(cmd == "clear")
		return input(ANSI_CLEAR);

	let output;
	try {
		// Get output from the command (only the synchronous commands that didn't use `&`).
		// Note: 2nd arg = callback that is called when an asynchronous command finishes.
		output = await $CLI.exec(cmd, out => $xterm.writeln(out));
		// Add extra break line so there's room to see what's going on in the terminal
	} catch (error) {
		output = error;
	}

	// Let parent component know we're done (used in Exercises to refresh status).
	// If useful, we can also send a 2nd arg containing data to send back.
	// We do this to keep Terminal.svelte independent from the `config.js` file
	// so it can be reused in other applications that don't have it.
	dispatch("status", "execDone");

	// Band-aid: don't show bowtie2 thread warnings
	if(typeof output === "string")
		output = output.replaceAll("pthread_sigmask() is not supported: this is a no-op.\n", "");

	// Ask the user for the next input
	return input(output);
}


// =============================================================================
// xterm.js handlers
// =============================================================================

// Keyboard shortcuts
function handleShortcuts(key)
{
	// Ctrl + L = Clear terminal (also clear input in case user had written something on the line)
	if(key.domEvent.ctrlKey && key.domEvent.key == "l") {
		$xterm.write(ANSI_CLEAR);
		$xtermAddons.echo.setInput("");
	}

	// Ctrl + A = Beginning of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "a")
		$xtermAddons.echo.setCursor(0);

	// Ctrl + E = End of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "e")
		$xtermAddons.echo.setCursor(Infinity);
}

// Autocomplete for subcommands. For some reason, it doesn't seem to work with local-echo. When I did "samtools <TAB>",
// it would output "samtools samtools <TAB>" along with the subcommands :(
async function handleAutocomplete(data)
{
	// Parse user input (only care about tabbing to autocomplete)
	if(data != "\t")
		return;
	const input = $xtermAddons.echo._input;
	const prgm = input.split(" ")[0].trim();  // ls /samtools/examples/toy
	const args = input.split(" ").slice(1);   // ["/samtools/examples/toy"]
	let cacheAutocomplete = [];

	// Turn off local-echo so we can override it's default behavior
	$xtermAddons.echo.detach();

	// Autocomplete main commands
	// e.g. "<TAB>" --> "samtools   bedtools"
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
			$xtermAddons.echo.handleCursorInsert(" ");
	}
	// If only one autocomplete option, then autocomplete it!
	else if(cacheAutocomplete.length == 1) {
		const remainingFragment = cacheAutocomplete[0].slice(args.length == 0 ? prgm.length : args[0].length);
		const extraSpacing = (prgm == "ls" && args.length != 0) ? "" : " ";
		$xtermAddons.echo.handleCursorInsert(remainingFragment + extraSpacing);
	// Otherwise, output all candidates
	} else {
		$xtermAddons.echo.printAndRestartPrompt(() => $xtermAddons.echo.printWide(cacheAutocomplete));
	}

	// Re-enable local-echo
	$xtermAddons.echo.attach();
}
</script>

<div bind:this={divTerminal} use:watchResize={handleResize} style="opacity: { ready ? 1 : 0.6 }; height:85vh; max-height:85vh; overflow:hidden">
	{#if !ready}
		<Spinner color="light" type="border" />
	{/if}
</div>
