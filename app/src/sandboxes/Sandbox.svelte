<script>
import { onMount } from 'svelte';

import { EditorView, basicSetup } from 'codemirror'
import { EditorState } from '@codemirror/state'

import { json } from "@codemirror/lang-json";
import { cpp } from "@codemirror/lang-cpp";
let elCodeMirror;
let editor;

onMount(() => {
	const data = {"abc": "sdf", "def": "fsd", "ghi": { "a": 123, "b": 456 }};
	initEditor(json, JSON.stringify(data, null, 2));
	console.log(editor);
});

function extract() {
	const code = editor.state.doc.toString();
	console.log(code);
}

function toggle() {
	initEditor(cpp, `function abc(){  return 123; }`);
}

function initEditor(lang, code) {
	// Define editor state
	const state = EditorState.create({
		doc: code,
		extensions: [
			basicSetup.slice(2),  // https://github.com/codemirror/basic-setup/blob/main/src/codemirror.ts#L47
			lang()
		]
	});

	// Creating editor for the first time
	if(!editor)
		editor = new EditorView({ state, parent: elCodeMirror });
	// Or update existing editor
	else
		editor.setState(state);
}
</script>

<div bind:this={elCodeMirror} style="border:1px solid #ccc" />

<button on:click={extract}>extract</button>

<button on:click={toggle}>toggle</button>
