<script>
import { onMount } from "svelte";
import { Terminal } from "xterm";
import { WebLinksAddon } from "xterm-addon-web-links";
import { SearchAddon } from "xterm-addon-search";
import LocalEchoController from "local-echo";
import "xterm/css/xterm.css";

// State
let term;
let termEcho;
let termSearch;
let divTerminal;

// -----------------------------------------------------------------------------
// On load
// -----------------------------------------------------------------------------
onMount(async () => {
	// Initialize terminal
	term = new Terminal()
	term.open(divTerminal);

	// Attach addons
	term.loadAddon(termEcho = new LocalEchoController());
	term.loadAddon(termSearch = new SearchAddon());
	term.loadAddon(new WebLinksAddon());

	// Register autocomplete handlers
	termEcho.addAutocompleteHandler(autocompleteCommands);

	// Welcome message
	term.writeln("# sandbox.bio\n");

	// Ready for user input
	term.focus();
	input();
})

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

// Auto-completes common commands
function autocompleteCommands(index, tokens)
{
	const command = tokens[0];
	console.log(index, tokens);

	// Root autocomplete
	if(index == 0)
		return ["samtools", "bedtools2"];
	// Samtools autocomplete. Need `&& tokens[index]`, otherwise results in "samtools samtools"
	if(index == 1 && command == "samtools" && tokens[index]) {
		return ["view", "index", "sort"];
	}
	return [];
}


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

async function exec(cmd)
{
	console.log(cmd);

	// 
	if(cmd == "clear")
		term.write("\x1bc");

	console.log(termSearch.findNext('sandbo'));
}

// Get user input
async function input()
{
	termEcho.read("$ ")
		.then(exec)
		.catch(error => console.error(`Error reading: ${error}`))
		.finally(input);
}
</script>

df
<div bind:this={divTerminal}></div>
