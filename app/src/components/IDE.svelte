<script>
export let input = "";

let divEditor;
let editor;
let pyodide;
let loading = {
	editor: false,
	pyodide: false,
};
let output = "";  // stdout, stderr
let result = "";  // return value of answer() function

// -----------------------------------------------------------------------------
// Run code
// -----------------------------------------------------------------------------

// Execute code with given input
function run() {
	output = "";
	try {
		pyodide.runPython(`${editor.getValue()}\n\nresult = dna_to_rna("${input}")`);
	} catch (error) {
		output += error;
	}
	result = pyodide.globals.get("result");
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
				value: `def dna_to_rna(t):\n\tprint('Output a string')\n\treturn t.replace('T', 'U')\n`,
				theme: "vs-light",
				language: "python",
				minimap: { enabled: false },
				automaticLayout: true
			});

			// Custom keyboard shortcuts
			editor.addAction({
				id: "execute-python",
				label: "Execute my script",
				keybindings: [ monaco.KeyMod.CtrlCmd | monaco.KeyCode.Enter ],
				// A precondition for this action.
				precondition: null,
				// A rule to evaluate on top of the precondition in order to dispatch the keybindings.
				keybindingContext: null,
				contextMenuGroupId: "navigation",
				contextMenuOrder: 1.5,
				// Method that will be executed when the action is triggered.
				run: run
			});
		});
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

<div>
	<div bind:this={divEditor} id="container-editor" class="border rounded-3 pt-2" style="height:50vh"></div>

	<div class="border rounded-3 p-2 mt-2" style="height:35vh; z-index:999; overflow-y:scroll">
		<h6>Input</h6>

		<div class="input-group mb-3">
			<input type="text" class="form-control font-monospace" id="input" bind:value={input}>
			<button class="btn btn-outline-secondary" type="button" on:click={run}>Run</button>
		</div>

		<h6 class="mt-3">
			Output
			<small class="text-muted" style="font-size:0.6em">Powered by <a href="https://pyodide.org/" target="_blank">Pyodide</a></small>
		</h6>
		<textarea id="result" class="form-control font-monospace" disabled>{result}</textarea>

		<h6 class="mt-3">Logs</h6>
		<textarea id="output" class="form-control font-monospace" disabled>{output}</textarea>
	</div>
</div>
