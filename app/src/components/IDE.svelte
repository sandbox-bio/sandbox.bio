<script>
import { createEventDispatcher } from "svelte";
import { EditorView, basicSetup } from "codemirror"
import { EditorState } from "@codemirror/state";
import { json } from "@codemirror/lang-json";
import { cpp } from "@codemirror/lang-cpp";

export let lang;  // json, cpp, null=no syntax highlighting
export let code;
export let editable = true;

// State
const dispatch = createEventDispatcher();
let divIDE;
let editor;

// (Re)initialize editor if the language changes (but not the code)
$: if(divIDE) initEditor(lang);

// If code changes, update the editor
$: if(editor && code !== editor.state.doc.toString()) {
	editor.dispatch({
		changes: {
			from: 0,
			to: editor.state.doc.length,
			insert: code
		}
	});
}

function initEditor(lang) {
	// Define basic extensions
	const extensions = [
		// Don't need all of basicSetup (https://github.com/codemirror/basic-setup/blob/main/src/codemirror.ts#L47)
		basicSetup.slice(2),
		// Support read-only editors if needed
		EditorView.editable.of(editable),
		// When the code is updated, ping parent component to update its state
		EditorView.updateListener.of(v => {
			if (v.docChanged)
				dispatch("update", v.state.doc.toString());
		})
	];

	// Add syntax highlighting
	if(lang === "cpp")
		extensions.push(cpp());
	else if(lang === "json")
		extensions.push(json());

	// Create or update editor
	const state = EditorState.create({
		doc: code,
		extensions,
		readOnly: true
	});
	if(!editor)
		editor = new EditorView({ state, parent: divIDE });
	// Or update existing editor
	else
		editor.setState(state);
}
</script>

<div bind:this={divIDE} />
