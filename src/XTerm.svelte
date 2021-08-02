<script>
import { onMount, createEventDispatcher } from "svelte";
import { Spinner } from "sveltestrap";
import "xterm/css/xterm.css";

// Imports
import { xterm, xtermAddons } from "terminal/xterm";
// import { CoreUtils } from "terminal/coreutils";
import { CLI } from "terminal/cli";

// Constants
const ANSI_CLEAR = "\x1bc";                // Clear terminal
const dispatch = createEventDispatcher();  // Dispatch for sending "exec" messages to parent component

// Autocomplete subcommands
const AUTOCOMPLETE = {
	samtools: () => ["view", "sort", "index", "idxstats"],
	bedtools: () => ["intersect", "merge", "complement", "genomecov", "jaccard", "makewindows", "flank"],
	ls: async args => {
		const pathSearch = args[0];                                               // /samtools/examples/toy
		const pathBase = pathSearch.substring(0, pathSearch.lastIndexOf("/")+1);  // /samtools/examples/
		// const files = await CoreUtils.ls([pathBase], true);
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
	rmdir: () => []
};


// =============================================================================
// State
// =============================================================================

export let ready = false;  // Whether CLI is ready for user input
export let intro;          // Intro string to display on Terminal once ready (optional)
let divTerminal;           // HTML element where terminal will be drawn
let runtime = 0;           // How long the last command took to run
let runtimeShow = false;   // Whether to output the command runtime

$: if(ready) input();      // Ask for user input once ready


// =============================================================================
// Initialization
// =============================================================================

onMount(() => {
	// Register handlers
	$xterm.onKey(handleShortcuts);
	$xterm.onData(handleAutocomplete);

	if(intro)
		$xterm.writeln(intro);

	// Prepare UI but don't allow input yet
	$xterm.open(divTerminal);
	$xtermAddons.fit.fit();
});


// =============================================================================
// xterm.js 
// =============================================================================

// Get user input 
function input(toPrint)
{
	runtime = window.performance.now() - runtime;
	if(toPrint)
		$xterm.writeln(toPrint);
	if(runtimeShow)
		$xterm.writeln(`\nRuntime: ${runtime}ms`);
	$xterm.focus();
	$xtermAddons.echo.read("$ ")
		.then(exec)
		.catch(console.error);
}

// Execute command
async function exec(cmd)
{
	console.log("[XTerm]", cmd);
	try {
		// Get synchronous outputs from the command
		const output = await $CLI.exec(cmd, out => {
			// This function is called when an asynchronous command finishes
			$xterm.writeln(out);
			// $xtermAddons.echo.writeln(out);
		});

		// Ask the user for the next input
		return input(output);
	} catch (error) {
		return input(error);
	}


	runtime = window.performance.now();
	runtimeShow = false;

	// -------------------------------------------------------------------------
	// Input validation
	// -------------------------------------------------------------------------
	if(!cmd || !(typeof cmd === "string"))
		return input();
	// Trim whitespace and remove duplicate spaces
	cmd = cmd.trim();
	cmd = cmd.replace(/  +/g, " ");
	cmd = cmd.replace(/\\/g, "");  // if multi-line command, remove \

	// Handle terminal-related commands
	if(cmd == "clear")
		return input(ANSI_CLEAR);

	// If user wants to check how long a command takes to run
	if(cmd.split(" ")[0] == "time") {
		runtimeShow = true;
		cmd = cmd.split(" ").slice(1).join(" ");
	}

	// -------------------------------------------------------------------------
	// Support basic file redirection (cmd > file) and piping (cmd | cmd2 | cmd3)
	// -------------------------------------------------------------------------

	// Prepare defaults
	let options = {
		cmd: cmd,
		file: null,  // null == stdout, otherwise save to a file path
		next: ""     // next command to run after this one (used by pipes)
	};

	// Redirections: assume ">" is not used in string arguments
	// e.g. samtools view -q20 toy.sam > toy.filtered.sam
	const redirections = cmd.split(">").map(d => d.trim());
	if(redirections.length > 2)
		return input("Unsupported command: Only support one '>' redirection.\n");
	else if(redirections.length == 2) {
		options.cmd = redirections[0];   // e.g. samtools view -q20 toy.sam
		options.file = redirections[1];  // e.g. toy.filtered.sam
	}

	// Piping: execute first part of the pipe, save output to a tmp file, then
	// run the other commands subsequently on those tmp files recursively.
	// e.g. samtools view | head -n 5 | wc -l > somefile
	const pipes = cmd.split("|").map(d => d.trim());
	if(pipes.length > 1)
	{
		// If we also have a redirection on top of the pipes, track it
		const fileRedir = options.file;  // somefile
		// And remove it from the command
		if(fileRedir)
			pipes[pipes.length-1] = pipes[pipes.length-1].split(">")[0];

		// Execute first part of the pipe, save the output to a temp file
		const fileTmp = `/tmp/tmp${parseInt(Math.random() * 10000)}`;
		options.cmd = pipes.shift();     // samtools view
		options.file = fileTmp;          // /tmp/tmp123

		// Next command should use the temp file we are writing to
		pipes[0] += ` ${fileTmp}`;       // head -n 5 /tmp/pipe123

		// Save remaining commands for the next time around
		options.next = pipes.join("|");  // head -n 5 /tmp/pipe123 | wc -l
		if(fileRedir)
			options.next += ` > ${fileRedir}`;  // head -n 5 /tmp/pipe123 | wc -l > somefile
	}

	// -------------------------------------------------------------------------
	// Execute command
	// -------------------------------------------------------------------------
	// Send message to parent component asking to execute the command.
	// The callback will output the result and ask for the next input.
	dispatch("exec", {
		cmd: options.cmd,
		callback: async out => {
			let stdout = "" + out;          // Convert to string if it's not already

			// Do we want to save this to a file?
			if(options.file) {
				// await CoreUtils.FS.writeFile(options.file, stdout);
				stdout = "";
			// Or just to stdout: append \n if needed
			} else if(!stdout.endsWith("\n")) {
				stdout += "\n";
			}

			// If we need to run more commands before we output, do so
			if(options.next)
				return await exec(options.next);
			// Otherwise, output to stdout and get next input
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

<div bind:this={divTerminal} style="opacity: { ready ? 1 : 0.6 }; height:85vh; max-height:85vh; overflow:hidden">
	{#if !ready}
		<Spinner color="light" type="border" />
	{/if}
</div>
