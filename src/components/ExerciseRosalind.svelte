<script>
import localforage from "localforage";
import { getLocalForageKey } from "$stores/config";
import { onMount } from "svelte";

export let expectedInput = ""; // Default input to show
export let expectedOutput = ""; // Expected output given that input
export let input = ""; // Current user input (starts out as expectedInput but can be modified by user)
export let code = ""; // Code to initialize the IDE with
export let fn = ""; // Function name; used to sync IDE state so use a unique name for this
export let fnParams = []; // List of function parameter names

// State
let divEditor;
let editor;
let pyodide;
let loaded = { editor: false, pyodide: false };
let output = ""; // stdout, stderr
let result = ""; // return value of answer() function
let success = null; // whether to show green or red border around output box (or nothing if null)
let updating = false; // if editor is being updated and we shouldn't save state
const CODE_LOADING = "Loading...";

// -----------------------------------------------------------------------------
// Reactive statements
// -----------------------------------------------------------------------------

// Initialize on mount
$: ready = loaded.editor === true && loaded.pyodide === true;
$: if (ready) editor.updateOptions({ readOnly: false });
$: {
	input = expectedInput;
	success = null;
}

// When code prop changes, what to do with existing code?
$: if (code && ready) updateEditor(code);

// -----------------------------------------------------------------------------
// Cache IDE state
// -----------------------------------------------------------------------------

async function updateEditor(newCode) {
	updating = true;

	// Check if there's something in localforage already?
	const data = await localforage.getItem(getLocalForageKey("ide") + fn);
	const codeIsDifferent = data !== null && data !== newCode;
	if (data !== null) newCode = data;

	// If not, update the editor
	editor.getModel().setValue(newCode);

	// If user made changes to default code, then run it (i.e. will show the correct success/fail colors)
	output = "";
	result = "";
	if (codeIsDifferent) run();

	updating = false;
}

// Save IDE state every few seconds
async function saveIDE(once = false) {
	if (ready && !updating) {
		console.log("Saving IDE state...");
		try {
			const state = editor.getValue();
			if (state != CODE_LOADING) await localforage.setItem(getLocalForageKey("ide") + fn, state);
		} catch (error) {}
	}
	if (once === false) setTimeout(saveIDE, 1500);
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
		if (fnParams.length > 1) {
			const paramsBreak = input.split("\n");
			const paramsSpaces = inputSanitized.split(" ");
			if (paramsBreak.length === fnParams.length) params = paramsBreak;
			else if (paramsSpaces.length === fnParams.length) params = paramsSpaces;
		}

		// Turn params into string or numbers
		params = params
			.map((d) => {
				if (isNumber(d)) return `${parseFloat(d)}`;
				return `"${d}"`;
			})
			.join(", ");

		// Run code
		pyodide.runPython(`${editor.getValue()}\n\nresult = ${fn}(${params})`);
	} catch (error) {
		output += error;
	}
	result = pyodide.globals.get("result");

	// Update success status
	if (input !== expectedInput || expectedOutput == "") success = null;
	else {
		let outputMatches = true;

		// Check that output is as expected, while supporting floating point output
		const observed = String(result).trim().split(/\n| /);
		const expected = expectedOutput.trim().split(/\n| /);

		if (observed.length != expected.length) outputMatches = false;
		else {
			for (let i = 0; i < observed.length; i++) {
				let valExpected = isNumber(expected[i]) ? parseFloat(expected[i]) : expected[i];
				let valObserved = isNumber(observed[i]) ? parseFloat(observed[i]) : observed[i];

				if (!isNumber(expected[i]) || !isNumber(observed[i])) {
					if (valObserved != valExpected) outputMatches = false;
				} else {
					const diff = Math.abs(valExpected - valObserved);
					if (diff >= 0.001) outputMatches = false;
				}
			}
		}

		success = outputMatches;
	}
}

// Reset IDE to original code
function reset() {
	updating = true;
	editor.getModel().setValue(code);
	updating = false;
}

// Utility function
function isNumber(d) {
	return !isNaN(d) && !isNaN(parseFloat(d));
}

// -----------------------------------------------------------------------------
// Initialization of Python + IDE
// -----------------------------------------------------------------------------

// Initialize Pyodide
async function initPython() {
	console.log("Initialize Python...");
	pyodide = await loadPyodide({
		indexURL: "https://cdn.jsdelivr.net/pyodide/v0.18.0/full/",
		stdout: (text) => (output += `${text}\n`),
		stderr: (text) => (output += `${text}\n`)
	});
	loaded.pyodide = true;
}

// Initialize Monaco Editor
async function initEditor() {
	console.log("Initialize editor...");

	// Note: "require" is provided by loader.min.js
	require.config({ paths: { vs: "https://cdnjs.cloudflare.com/ajax/libs/monaco-editor/0.26.1/min/vs" } });
	require(["vs/editor/editor.main"], () => {
		editor = monaco.editor.create(divEditor, {
			value: CODE_LOADING,
			theme: "vs-light",
			language: "python",
			minimap: { enabled: false },
			automaticLayout: true,
			readOnly: true
		});

		// Custom IDE menu items
		editor.addAction({
			id: "execute-python",
			label: "Execute my script",
			keybindings: [monaco.KeyMod.CtrlCmd | monaco.KeyCode.Enter],
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

		loaded.editor = true;
	});
}

// Import 3rd-party tools in the right order
onMount(async () => {
	function loadScript(url) {
		return new Promise(function (resolve, reject) {
			let script = document.createElement("script");
			script.src = url;
			script.async = false;
			script.onload = () => resolve(url);
			script.onerror = () => reject(url);
			document.body.appendChild(script);
		});
	}

	// The order matters here for some reason...
	await loadScript("https://cdn.jsdelivr.net/pyodide/v0.18.0/full/pyodide.js").then(initPython);
	await loadScript("https://cdnjs.cloudflare.com/ajax/libs/monaco-editor/0.26.1/min/vs/loader.min.js").then(initEditor);
});
</script>

<svelte:head>
	<!-- For some reason, loading both monaco and pyodide here fails when monaco loads first... -->
	<!-- <script async src="https://cdn.jsdelivr.net/pyodide/v0.18.0/full/pyodide.js" on:load={initPython}></script> -->
	<!-- <script src="https://cdnjs.cloudflare.com/ajax/libs/monaco-editor/0.26.1/min/vs/loader.min.js" on:load={initEditor}></script> -->
</svelte:head>

<div>
	<div bind:this={divEditor} class="border rounded-3 pt-2" style="height:50vh" />

	<div class="border rounded-3 p-2 mt-2" style="height:35vh; z-index:999; overflow-y:scroll">
		<h6>Input</h6>

		<div class="input-group mb-3">
			<textarea
				id="sdf"
				class="form-control font-monospace"
				disabled={!ready}
				bind:value={input}
				on:keydown={(e) => {
					if (e.key === "Enter" && e.metaKey === true) run();
				}}
				rows="2"
			/>
			<button class="btn btn-primary" type="button" disabled={!ready} on:click={run}>Run</button>
		</div>

		<h6 class="mt-3">Output</h6>
		<textarea
			id="result"
			class="form-control font-monospace"
			class:border-2={success !== null}
			class:border-success={success === true}
			class:border-danger={success === false}
			rows="2"
			disabled>{result}</textarea
		>

		<h6 class="mt-3">Logs</h6>
		<textarea id="output" class="form-control font-monospace" rows="2" disabled>{output}</textarea>
	</div>
</div>

<style>
textarea {
	font-size: 0.75em;
}
</style>
