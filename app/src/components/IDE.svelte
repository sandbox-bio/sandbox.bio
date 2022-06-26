<script>
import { onMount } from "svelte";
import { createEventDispatcher } from "svelte";
import { EditorView, basicSetup } from "codemirror"
import { EditorState } from "@codemirror/state";
import { json } from "@codemirror/lang-json";  // json input/output
import { cpp } from "@codemirror/lang-cpp";  // cmd

export let lang = null;  // json, cpp, null=no syntax highlighting
export let code = `"sup"`;

const dispatch = createEventDispatcher();
let divIDE;
let editor;

$: if(divIDE && lang) {
	console.log(editor);
	initEditor(lang);
}

// onMount(() => {
// 	initEditor(lang, code);
// 	// const data = {"abc": "sdf", "def": "fsd", "ghi": { "a": 123, "b": 456 }};
// 	// initEditor(json, JSON.stringify(data, null, 2));
// 	// console.log(editor);
// });

// function extract() {
// 	const code = editor.state.doc.toString();
// 	console.log(code);
// }

// function toggle() {
// 	initEditor(cpp, `function abc(){  return 123; }`);
// }

function initEditor(lang, doc) {
	// Define editor state
	const extensions = [
		// Don't need all of basicSetup (https://github.com/codemirror/basic-setup/blob/main/src/codemirror.ts#L47)
		basicSetup.slice(2),
		// When the code is updated, ping parent component to update its state
		EditorView.updateListener.of(v => {
			if (v.docChanged)
				dispatch("update", v.state.doc.toString());
		})
	];  
	if(lang === "cpp")
		extensions.push(cpp());
	else if(lang === "json")
		extensions.push(json());
	const state = EditorState.create({
		doc: doc || code,
		extensions
	});

	// Creating editor for the first time
	if(!editor)
		editor = new EditorView({ state, parent: divIDE });
	// Or update existing editor
	else
		editor.setState(state);
}
</script>

<div bind:this={divIDE} />
