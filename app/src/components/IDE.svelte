<script>
import localforage from "localforage";
import { getLocalForageKey } from "stores/config";

export let expectedInput = "";   // Default input to show
export let expectedOutput = "";  // Expected output given that input
export let input = "";           // Current user input (starts out as expectedInput but can be modified by user)
export let code = "";            // Code to initialize the IDE with
export let fn = "";              // Function name; used to sync IDE state so use a unique name for this
export let fnParams = [];        // List of function parameter names

// State
let divEditor;
let editor;
let pyodide;
let loading = { editor: false, pyodide: false };
let loaded = { editor: false, pyodide: false };
let output = "";       // stdout, stderr
let result = "";       // return value of answer() function
let success = null;    // whether to show green or red border around output box (or nothing if null)
let updating = false;  // if editor is being updated and we shouldn't save state
const CODE_LOADING = "Loading...";


// -----------------------------------------------------------------------------
// Reactive statements
// -----------------------------------------------------------------------------

// Initialize on mount
$: ready = loaded.editor === true && loaded.pyodide === true;
$: if(ready) editor.updateOptions({ readOnly: false });
$: { input = expectedInput; success = null; };

// When code prop changes, what to do with existing code?
$: if(code && ready) updateEditor(code);


// -----------------------------------------------------------------------------
// Cache IDE state
// -----------------------------------------------------------------------------

async function updateEditor(newCode) {
	updating = true;

	// Check if there's something in localforage already?
	const data = await localforage.getItem(getLocalForageKey("ide") + fn);
	const codeIsDifferent = data !== null && data !== newCode;
	if(data !== null)
		newCode = data;

	// If not, update the editor
	editor.getModel().setValue(newCode);

	// If user made changes to default code, then run it (i.e. will show the correct success/fail colors)
	output = "";
	result = "";
	if(codeIsDifferent)
		run();

	updating = false;
}

// Save IDE state every few seconds
async function saveIDE(once=false) {
	if(ready && !updating) {
		console.log("Saving IDE state...");
		try {
			const state = editor.getValue();
			if(state != CODE_LOADING)
				await localforage.setItem(getLocalForageKey("ide") + fn, state);
		} catch (error) {}
	}
	if(once === false)
		setTimeout(saveIDE, 1500);
}
saveIDE();


// -----------------------------------------------------------------------------
// Run code
// -----------------------------------------------------------------------------

// Execute code with given input
function run() {
	output = "";
	try {
		// Figure out function parameters
		const inputSanitized = input.replaceAll("\n", "\\n");
		// By default, just pass one argument
		let params = [inputSanitized];
		// But could need more than 1
		if(fnParams.length > 1) {
			const paramsBreak = input.split("\n");
			const paramsSpaces = inputSanitized.split(" ");
			if(paramsBreak.length === fnParams.length)
				params = paramsBreak;
			else if(paramsSpaces.length === fnParams.length)
				params = paramsSpaces;
		}

		// Turn params into string or numbers
		params = params.map(d => {
			if(isNumber(d))
				return `${parseFloat(d)}`;
			return `"${d}"`;
		}).join(", ");

		// Run code
		pyodide.runPython(`${editor.getValue()}\n\nresult = ${fn}(${params})`);
	} catch (error) {
		output += error;
	}
	result = pyodide.globals.get("result");

	// Update success status
	if(input !== expectedInput || expectedOutput == "")
		success = null;
	else {
		let outputMatches = true;

		// Check that output is as expected, while supporting floating point output
		const observed = String(result).trim().split(/\n| /);
		const expected = expectedOutput.trim().split(/\n| /);

		if(observed.length != expected.length)
			outputMatches = false;
		else {
			for(let i = 0; i < observed.length; i++) {
				let valExpected = isNumber(expected[i]) ? parseFloat(expected[i]) : expected[i];
				let valObserved = isNumber(observed[i]) ? parseFloat(observed[i]) : observed[i];

				if(!isNumber(expected[i]) || !isNumber(observed[i])) {
					if(valObserved != valExpected)	
						outputMatches = false;
				} else {
					const diff = Math.abs(valExpected - valObserved);
					if(diff >= 0.001)
						outputMatches = false;
				}
			}
		}

		success = outputMatches;
	}
}

function isNumber(d) {
	return !isNaN(d) && !isNaN(parseFloat(d));
}

// Reset IDE to original code
function reset() {
	updating = true;
	editor.getModel().setValue(code);
	updating = false;
}


// -----------------------------------------------------------------------------
// Initialization of Python + IDE
// -----------------------------------------------------------------------------

// Initialize Pyodide
async function initPython(){
	console.log("Initialize Python...");
	if(loading.pyodide)
		return;
	loading.pyodide = true;

	try {
		pyodide = await loadPyodide({
			indexURL : "https://cdn.jsdelivr.net/pyodide/v0.18.0/full/",
			stdout: text => output += `${text}\n`,
			stderr: text => output += `${text}\n`,
		});
		loaded.pyodide = true;
	} catch (error) {
		console.log("Pyodide Failed");
		loading.pyodide = false;
	}
}

// Initialize Monaco Editor
async function initEditor()
{
	console.log("Initialize editor...");
	if(loading.editor)
		return;
	loading.editor = true;

	try {
		// require is provided by loader.min.js
		require.config({ paths: { vs: "https://cdnjs.cloudflare.com/ajax/libs/monaco-editor/0.26.1/min/vs" }});
		require(["vs/editor/editor.main"], () => {
			editor = monaco.editor.create(divEditor, {
				value: CODE_LOADING,
				theme: "vs-light",
				language: "python",
				minimap: { enabled: false },
				automaticLayout: true,
				readOnly: true
			});

			// Custom keyboard shortcuts
			editor.addAction({
				id: "execute-python",
				label: "Execute my script",
				keybindings: [ monaco.KeyMod.CtrlCmd | monaco.KeyCode.Enter ],
				contextMenuGroupId: "navigation",
				contextMenuOrder: 1.5,
				run: run
			});

			editor.addAction({
				id: "reset-script",
				label: "Reset code to default (changes will be lost)",
				contextMenuGroupId: "navigation",
				contextMenuOrder: 1.5,
				run: reset
			});
		});
		loaded.editor = true;
	} catch (error) {
		console.log("Editor Failed");
		loading.editor = false;
	}
}

// Initialize on page load
async function init() {
	console.log("init")
	if(!loading.pyodide)
		initPython();
	if(!loading.editor)
		initEditor();
	if(!loading.pyodide || !loading.editor)
		setTimeout(init, 400);
}
init();
</script>

<svelte:head>
	<script src="https://cdn.jsdelivr.net/pyodide/v0.18.0/full/pyodide.js"></script>
	<script src="https://cdnjs.cloudflare.com/ajax/libs/monaco-editor/0.26.1/min/vs/loader.min.js"></script>
</svelte:head>

<style>
textarea {
	font-size: 0.75em;
}
</style>

<div>
	<div bind:this={divEditor} id="container-editor" class="border rounded-3 pt-2" style="height:50vh"></div>

	<div class="border rounded-3 p-2 mt-2" style="height:35vh; z-index:999; overflow-y:scroll">
		<h6>Input</h6>

		<div class="input-group mb-3">
			<textarea id="sdf" class="form-control font-monospace" disabled={!ready} bind:value={input} on:keydown={e => {
				if(e.key === "Enter" && e.metaKey === true)
					run();
			}} rows="2"></textarea>
			<button class="btn btn-primary" type="button" disabled={!ready} on:click={run}>Run</button>
		</div>

		<h6 class="mt-3">
			Output
			<small class="text-muted" style="font-size:0.6em">Powered by <a href="https://pyodide.org/" target="_blank">Pyodide</a></small>
		</h6>
		<textarea id="result" class="form-control font-monospace" class:border-2={success !== null} class:border-success={success === true} class:border-danger={success === false} rows="2" disabled>{result}</textarea>

		<h6 class="mt-3">Logs</h6>
		<textarea id="output" class="form-control font-monospace" rows="2" disabled>{output}</textarea>
	</div>
</div>
