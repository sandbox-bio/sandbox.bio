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
const TOOLS_DEFAULT = ["samtools/1.10", "bcftools/1.10", "bedtools/2.29.2", "bowtie2/bowtie2-align-s/2.4.2"];

// Autocomplete subcommands
const AUTOCOMPLETE = {
	samtools: ["view", "sort", "depth", "index", "idxstats", "flags", "flagstats"],
	bedtools: ["intersect", "merge", "complement", "genomecov", "jaccard", "makewindows", "flank"],
	bcftools: ["view", "index", "call", "query", "merge"],
	bowtie2: [],
	ls: [],
	ll: [],
	cat: [],
	head: [],
	tail: [],
	wc: [],
	pwd: [],
	cd: [],
	echo: [],
	mv: [],
	rm: [],
	mkdir: [],
	rmdir: [],
	env: []
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

	// Delete autocomplete programs if we haven't loaded them
	try {
		TOOLS_DEFAULT.map(toolName => {
			if(tools.find(d => d == toolName) == null)
				delete AUTOCOMPLETE[toolName.split("/")[0]];
		});
	} catch (error) {}

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
		output = await $CLI.exec(cmd, out => $xterm.write(out));
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
		const originalInput = $xtermAddons.echo._input;
		$xterm.write(ANSI_CLEAR);
		$xtermAddons.echo.setInput(originalInput);
	}

	// Ctrl + A = Beginning of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "a")
		$xtermAddons.echo.setCursor(0);

	// Ctrl + E = End of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "e")
		$xtermAddons.echo.setCursor(Infinity);

	// Ctrl + W = Delete last word
	if(key.domEvent.ctrlKey && key.domEvent.key == "w")
		$xtermAddons.echo.handleData("\x1b\x7F");
}

// Autocomplete for subcommands. For some reason, it doesn't seem to work with local-echo. When I did "samtools <TAB>",
// it would output "samtools samtools <TAB>" along with the subcommands :(
// NOTE: this only supports autocompleting the ends of commands, not in the middle
async function handleAutocomplete(data)
{
	// Parse user input (only care about tabbing to autocomplete)
	if(data != "\t")
		return;
	const input = $xtermAddons.echo._input;
	const prgm = input.split(" ")[0].trim();      // ls /samtools/examples/toy
	const args = input.split(" ").slice(1);       // ["/samtools/examples/toy"]
	const userFragment = input.split(" ").pop();  // The fragment the user entered
	let cacheAutocomplete = [];                   // Options to list to the user
	let appendExtraSpace = true;                  // Whether we should append a space to the user's input (do that if have 1 option left and it's not a folder listing)

	// Turn off local-echo so we can override it's default behavior
	$xtermAddons.echo.detach();

	// Autocomplete main commands
	// e.g. "<TAB>" --> "samtools   bedtools"
	if(args.length == 0)
		cacheAutocomplete = Object.keys(AUTOCOMPLETE);

	// Autocomplete subcommands
	else if(prgm in AUTOCOMPLETE)
	{
		// e.g. "samtools " --> "view   sort   index  idxstats"
		// e.g. "samtools i" --> "index   idxstats"
		if(args.length < 2)
			cacheAutocomplete = AUTOCOMPLETE[prgm];

		// Autocomplete variables and file paths if no subcommands available
		// e.g. "samtools view to" --> "samtools view toy.bam"
		if(cacheAutocomplete.length == 0) {
			// Autocomplete variables
			if(userFragment.startsWith("$")) {
				cacheAutocomplete = Object.keys($CLI._vars).map(d => `$${d}`);

			// Autocomplete file paths
			} else {
				// Infer base path and files within it (default to `.`)
				const pathBase = userFragment.substring(0, userFragment.lastIndexOf("/") + 1);
				const files = await $CLI.coreutils.ls([ pathBase || "." ], true);
				// Prepend base path since `ls` doesn't do that for us
				cacheAutocomplete = files.map(d => pathBase + d.name);
			}
		}
	}

	// Only show those that match what the user entered
	cacheAutocomplete = cacheAutocomplete.filter(d => d.startsWith(userFragment));
	// If we're listing folders and there's only 1 option, don't add a space yet
	// because the user will want to keep autocompleting within that folder!
	if(cacheAutocomplete.length == 1 && cacheAutocomplete[0].endsWith("/"))
		appendExtraSpace = false;

	// Process autocomplete
	// If nothing to autocomplete, just add a space
	if(cacheAutocomplete.length == 0) {
		if(!input.endsWith(" "))
			$xtermAddons.echo.handleCursorInsert(" ");

	// If only one autocomplete option, then autocomplete it!
	} else if(cacheAutocomplete.length == 1) {
		const remainingFragment = cacheAutocomplete[0].slice(userFragment.length);
		$xtermAddons.echo.handleCursorInsert(remainingFragment + (appendExtraSpace ? " " : ""));

	// Otherwise, output all candidates
	} else {
		// Do all remaining candidates share a substring in common? If so, output it
		const sharedFragment = getSharedSubstring(cacheAutocomplete);
		$xtermAddons.echo.handleCursorInsert(sharedFragment.replace(userFragment, ""));
		// Output all candidates
		$xtermAddons.echo.printAndRestartPrompt(() => $xtermAddons.echo.printWide(cacheAutocomplete));
	}

	// Re-enable local-echo
	$xtermAddons.echo.attach();
}

// Source: <https://stackoverflow.com/a/1917041>
function getSharedSubstring(array){
	const A = array.concat().sort();
	let a1 = A[0], a2 = A[A.length-1], L = a1.length, i = 0;
	while(i < L && a1.charAt(i) === a2.charAt(i)) i++;
	return a1.substring(0, i);
}
</script>

<div bind:this={divTerminal} use:watchResize={handleResize} style="opacity: { ready ? 1 : 0.6 }; height:85vh; max-height:85vh; overflow:hidden">
	{#if !ready}
		<Spinner color="light" type="border" />
	{/if}
</div>
