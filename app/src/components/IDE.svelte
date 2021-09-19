<script>
import { TabContent, TabPane } from "sveltestrap";

export let input = "";

let divEditor;
let editor;
let pyodide;
let loading = {
	editor: false,
	pyodide: false,
};
let output = "";

async function initPython(){
	console.log("Initialize Python...")
	if(loading.pyodide)
		return;
	loading.pyodide = true;

	try {
		pyodide = await loadPyodide({
			indexURL : "https://cdn.jsdelivr.net/pyodide/v0.18.0/full/",
			stdout: (text) => output += text,
			stderr: (text) => output += text,
		});
		console.log(pyodide.runPython("1 + 2"));
	} catch (error) {
		console.log("Pyodide Failed");
		loading.pyodide = false;
	}
}

async function initEditor()
{
	console.log("Initialize editor...")
	if(loading.editor)
		return;
	loading.editor = true;

	try {
		// require is provided by loader.min.js
		require.config({ paths: { vs: "https://cdnjs.cloudflare.com/ajax/libs/monaco-editor/0.26.1/min/vs" }});
		require(["vs/editor/editor.main"], () => {
			editor = monaco.editor.create(divEditor, {
				value: `def answer(t):\n\treturn t.replace('T', 'U')\n`,
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
				run: function(ed) {
					output = "";
					try {
						output = pyodide.runPython(ed.getValue() + `\nanswer("${input}")`);
					} catch (error) {
						output = error;
					}
					return null;
				}
			});
		});
	} catch (error) {
		console.log("Editor Failed");
		loading.editor = false;
	}
}

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

	<div bind:this={divEditor} id="container-editor" class="border rounded-3 pt-2" style="height:60vh"></div>

	<div class="border rounded-3 p-2 mt-2" style="height:25vh; z-index:999; overflow-y:scroll">
		<TabContent>
			<TabPane tabId="input" tab="Input" active>
				<textarea class="form-control m-2" bind:value={input}></textarea>
			</TabPane>
			<TabPane tabId="output" tab="Output">
				<pre class="m-2">{output}</pre>
			</TabPane>
		</TabContent>
		  
	

		<!-- <h5>
			Output
			<small class="text-muted" style="font-size:0.6em">Powered by <a href="https://pyodide.org/" target="_blank">PyIodide</a></small>
		</h5>
		<pre>{output}</pre> -->
	</div>
</div>
