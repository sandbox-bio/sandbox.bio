<script>
import { onMount, createEventDispatcher } from "svelte";
import { watchResize } from "svelte-watch-resize";
import { Spinner, Table, Modal } from "sveltestrap";
import AnsiUp from "ansi_up";
import "xterm/css/xterm.css";

// Imports
import { xterm, xtermAddons } from "terminal/xterm";
import { CLI } from "terminal/cli";
import { config, env, MAX_FILE_SIZE_TO_CACHE } from "./stores/config";
import { status } from "./stores/status";
import { tutorial } from "./stores/tutorial";

// Constants
const ANSI_CLEAR = "\x1bc";
const COREUTILS = ["basename", "cat", "comm", "cut", "date", "echo", "fold", "head", "join", "ls", "md5sum", "paste", "seq", "shuf", "sort", "tail", "tee", "tr", "uniq", "wc"];
const HTSLIB_TOOLS = ["tabix", "htsfile", "bgzip"];
const TOOLS_DEFAULT = [
	// Bioinformatics
	{ loading: "lazy", tool: "samtools", version: "1.10" },
	{ loading: "lazy", tool: "bcftools", program: "bcftools", version: "1.10" },
	{ loading: "lazy", tool: "bedtools", version: "2.29.2" },
	{ loading: "lazy", tool: "bowtie2", program: "bowtie2-align-s", version: "2.4.2" },
	{ loading: "lazy", tool: "minimap2", version: "2.22" },
	{ loading: "lazy", tool: "seqtk", version: "1.3" },
	{ loading: "lazy", tool: "kalign", version: "3.3.1" },
	{ loading: "lazy", tool: "ivar", version: "1.3.1" },
	{ loading: "lazy", tool: "fasttree", version: "2.1.11" },
	{ loading: "lazy", tool: "fastp", version: "0.20.1" },
	{ loading: "lazy", tool: "tree", version: "2.0.4" },
	{ loading: "lazy", tool: "lsd2", version: "2.3" },
	{ loading: "lazy", tool: "tn93", version: "1.0.9" },
	...HTSLIB_TOOLS.map(program => ({ program, tool: "htslib", loading: "lazy", version: "1.10", reinit: true })),

	// General terminal tools
	{ loading: "lazy", tool: "jq", version: "1.6" },
	{ loading: "lazy", tool: "gawk", version: "5.1.0", reinit: true },
	{ loading: "lazy", tool: "grep", version: "3.7", reinit: true },
	{ loading: "lazy", tool: "sed", version: "4.8", reinit: true },
	{ loading: "lazy", tool: "findutils", program: "find", version: "4.9.0", reinit: true },
	...COREUTILS.map(program => ({ program, tool: "coreutils", loading: "lazy", version: "8.32", reinit: true }))
];

// Generate default autocomplete
const AUTOCOMPLETE_DEFAULT = {};
TOOLS_DEFAULT.forEach(d => {
	if(d.tool === "bowtie2")
		AUTOCOMPLETE_DEFAULT[d.tool] = [];
	else
		AUTOCOMPLETE_DEFAULT[d.program || d.tool] = [];
});

// Autocomplete subcommands
const AUTOCOMPLETE = {
	...AUTOCOMPLETE_DEFAULT,
	// Bioinformatics tools
	samtools: ["view", "sort", "depth", "index", "idxstats", "flags", "flagstats"],
	bedtools: ["intersect", "merge", "complement", "genomecov", "jaccard", "makewindows", "flank"],
	bcftools: ["view", "index", "call", "query", "merge"],
	ivar: ["trim", "variants", "filtervariants", "consensus", "getmasked", "removereads", "version"],
	seqtk: ["seq", "comp", "sample", "subseq", "fqchk", "mergepe", "trimfq", "hety", "gc", "mutfa", "mergefa", "famask", "dropse", "rename", "randbase", "cutN", "listhet"],
	// Open/download files
	open: [], download: [],
	// Host info
	hostname: [], uname: [], whoami: [], man: [],
	// Env variables
	env: [], unset: [], history: [],
	// Aliases
	ll: [],
	// Simulated coreutils
	pwd: [], cd: [], mv: [], rm: [], cp: [], mkdir: [], rmdir: [], touch: [],
};


// =============================================================================
// State
// =============================================================================

export let intro = "";                     // Intro string to display on Terminal once ready (optional)
export let init = "";                      // Command to run to initialize the environment (optional)
export let files = [];                     // Files to preload on the filesystem
export let tools = TOOLS_DEFAULT;          // Aioli tools to load
export let pwd = "";                       // Path relative to /shared/data where user should start at

let aioliReady = false;                    // Equals true once Aioli is initialized
let divTerminal;                           // HTML element where terminal will be drawn
let fileInput;                             // Hidden HTML file input element for mounting local file
let nbInit = 0;                            // Number of times we've reinitialized the terminal (i.e. when user logs in/out)
let modalKbdOpen = false;                  // Set to true when the shortcuts modal is open
let modalKbdToggle = () => modalKbdOpen = !modalKbdOpen;
const dispatch = createEventDispatcher();  // Send info to parent component when cmd is done

$: ready = aioliReady && $status.app;      // Ready to go once both Aioli and the app are initialized
$: if(ready) initTerminal();               // Ask for user input once ready


// =============================================================================
// Initialization
// =============================================================================

// Load filesystem from cache and get user input
async function initTerminal() {
	await $CLI.fsLoad($tutorial);
	nbInit++;
	input();
	saveFS();
}

// Save filesystem state every few seconds
async function saveFS() {
	try {
		await $CLI.fsSave($tutorial);
		console.log("Saved FS state");
	} catch (error) {}
	setTimeout(saveFS, 1000);
}

// On mount
onMount(async () => {
	// Register handlers
	$xterm.onKey(handleShortcuts);
	$xterm.onData(handleAutocomplete);

	// Write out an intro if any specified
	if(intro)
		$xterm.writeln(intro);

	// Prepare UI but don't allow input yet
	$xterm.open(divTerminal);

	// Put main tools at the beginning and the rest at the end of the list of tools.
	// That way, we eager load the important stuff but lazy load the rest.
	const toolsEagerLoad = TOOLS_DEFAULT
		.filter(tool => tools.includes(tool.tool) || tools.includes(tool.program))
		.map(tool => { tool.loading = "eager"; return tool; });
	const toolsLazyLoad = TOOLS_DEFAULT.filter(tool => !tools.includes(tool.tool) && !tools.includes(tool.program));
	tools = [...toolsEagerLoad, ...toolsLazyLoad];

	// Initialize Aioli
	try {
		await $CLI.init({ tools, files, pwd });
		// Custom command to run once terminal is ready
		if(init)
			await $CLI.exec(init);
		aioliReady = true;
	} catch (error) {
		console.error("Could not load terminal");
		console.error(error)
	}
});

// Resize xterm when the window size changes
function handleResize() {
	$xtermAddons.fit.fit();
}

// =============================================================================
// Sidebar operations
// =============================================================================

// Export ANSI to HTML and open in new tab
function exportTerminal() {
	const terminalRaw = $xtermAddons.serialize.serialize().replace(/\[[0-9]C/g, "\t");
	const terminalHTML = "<pre>" + (new AnsiUp).ansi_to_html(terminalRaw) + "</pre>";
	const blob = new Blob([ terminalHTML ], { type: "text/html" });
	const url = URL.createObjectURL(blob);
	window.open(url);
}

// Mount local file to virtual file system
async function mountLocalFile(event) {
	const files = event.target.files;
	if(!files) {
		console.warn("No file specified.");
		return;
	}

	// Note that files that already exist will be overwritten!
	const paths = await $CLI.utils.mount(files);
	const pathsTxt = paths.map((path, i) => {
		if(!files[i]) return;
		const extra = files[i].size <= MAX_FILE_SIZE_TO_CACHE ? "" : "   (large file; won't persist on page refresh)";
		return `#   ${path}${extra}`;
	}).join("\n");
	input(`\n\u001b[0;32m# Files mounted:\n${pathsTxt}\u001b[0m\n\n`);
}

// Clear command-line history
async function clearHistory() {
	$xterm.writeln("history -c");
	exec("history -c");
}


// =============================================================================
// xterm.js 
// =============================================================================

// Get user input 
async function input(toPrint)
{
	if(nbInit > 1) {
		$xterm.writeln("\n");
		nbInit = 1;
	}
	if(toPrint)
		$xterm.write(toPrint);
	$xterm.focus();

	// Prepare prompt, e.g. "guest@sandbox$ "
	const prompt = $env["PS1"].replaceAll('\\u', $env["USER"]).replaceAll('\\h', $config.hostname);
	$xtermAddons.echo.read(`\u001b[1;34m${prompt}\u001b[0m`)
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
		output = await $CLI.exec(cmd, out => $xterm.write(filter(out)));
		// Add extra break line so there's room to see what's going on in the terminal
	} catch (error) {
		output = error;
	}

	// Analytics
	try {
		const isTutorial = $tutorial?.id && $tutorial.id !== "playground";
		fetch(`${$config.api}/ping`, {
			method: "POST",
			mode: "no-cors",
			body: JSON.stringify({
				playground: isTutorial ? null : "terminal",
				tutorial: isTutorial ? $tutorial.id : null,
				step_from: isTutorial ? $tutorial.step : null,
				step_to: isTutorial ? $tutorial.step : null,
			})
		});
	} catch (error) {}

	// Let parent component know we're done (used in Exercises to refresh status).
	// If useful, we can also send a 2nd arg containing data to send back.
	// We do this to keep Terminal.svelte independent from the `config.js` file
	// so it can be reused in other applications that don't have it.
	dispatch("status", "execDone");

	// Ask the user for the next input
	return input(output);
}

// Filter out warnings
function filter(output) {
	// Band-aid: don't show bowtie2 thread warnings
	if(typeof output === "string") {
		output = output.replaceAll("pthread_sigmask() is not supported: this is a no-op.\n", "");
		output = output.replaceAll("warning: unsupported syscall: __sys_prlimit64\n", "");
	}

	return output;
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
				cacheAutocomplete = Object.keys($env).map(d => `$${d}`);

			// Autocomplete file paths
			} else {
				// Infer base path and files within it (default to `.`)
				const pathBase = userFragment.substring(0, userFragment.lastIndexOf("/") + 1);
				const files = await $CLI.utils.ls([ pathBase || "." ], true);
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


<!-- Terminal -->
<div bind:this={divTerminal} use:watchResize={handleResize} style="opacity: { ready ? 1 : 0.6 }; height:85vh; max-height:85vh; overflow:hidden">
	<!-- Hamburger menu for settings -->
	<div class="cli-options text-muted">
		<button class="btn btn-outline-secondary p-0 m-0 border-0" type="button" data-bs-toggle="dropdown" aria-expanded="false">
			<i class="bi bi-three-dots-vertical"></i>
		</button>
		<ul class="dropdown-menu">
			<li><button class="dropdown-item" on:click={() => fileInput.click()}>Mount local files</button></li>
			<li><button class="dropdown-item" on:click={exportTerminal}>Export as HTML</button></li>
			<li><button class="dropdown-item" on:click={clearHistory}>Clear command history</button></li>
			<li><button class="dropdown-item" on:click={modalKbdToggle}>Keyboard Shortcuts</button></li>
		</ul>
	</div>
	{#if !ready}
		<Spinner color="light" type="border" />
	{/if}
</div>

<!-- Hidden input file for mounting local files -->
<input type="file" on:change={mountLocalFile} bind:this={fileInput} style="display:none" multiple />

<!-- Keyboard Shortcuts Modal -->
<Modal body header="Keyboard Shortcuts" isOpen={modalKbdOpen} toggle={modalKbdToggle}>
	<Table>
		<thead>
			<tr>
				<th>Shortcut</th>
				<th>Action</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td><code>Ctrl + L</code></td>
				<td>Clear</td>
			</tr>
			<tr>
				<td><code>Ctrl + A</code></td>
				<td>Go to start of line</td>
			</tr>
			<tr>
				<td><code>Ctrl + E</code></td>
				<td>Go to end of line</td>
			</tr>
			<tr>
				<td><code>Ctrl + W</code></td>
				<td>Delete previous word</td>
			</tr>
			<tr>
				<td><code>Alt + Left</code></td>
				<td>Go to previous word</td>
			</tr>
			<tr>
				<td><code>Alt + Right</code></td>
				<td>Go to following word</td>
			</tr>
		</tbody>
	</Table>
</Modal>

<style>
/* Hamburger menu */
.cli-options {
	position: absolute;
	right: 0;
	z-index: 100;
	cursor: pointer;
}

.cli-options:hover {
	color: white !important;
}
</style>
