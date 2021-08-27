<script>
// Had to modify rollup.config.js to set inlineDynamicImports to true:
// 		export default {
// 			input: "src/main.js",
// 			output: {
// 				sourcemap: true,
// 				format: "iife",
// 				name: "app",
// 				file: "public/build/bundle.js",
// 				inlineDynamicImports: true
// 			},
// ---------------------------------------------------------------

import { onMount } from "svelte";

// ---------------------------------------------------------------
// // Size of total build/bundle.js: 15MB
// ---------------------------------------------------------------
// import * as monaco from 'monaco-editor';

// ---------------------------------------------------------------
// Size of total build/bundle.js: 5.8MB, 2.2MB bundle
// ---------------------------------------------------------------
// import 'monaco-editor/esm/vs/base/browser/ui/';
// import 'monaco-editor/esm/vs/editor/browser/controller/coreCommands.js';
// import 'monaco-editor/esm/vs/editor/browser/widget/codeEditorWidget.js';
import 'monaco-editor/esm/vs/basic-languages/python/python.contribution.js';
// import 'monaco-editor/esm/vs/editor/contrib/bracketMatching/bracketMatching.js';
// import 'monaco-editor/esm/vs/editor/contrib/clipboard/clipboard.js';
// import 'monaco-editor/esm/vs/editor/contrib/comment/comment.js';
// import 'monaco-editor/esm/vs/editor/contrib/contextmenu/contextmenu.js';
// import 'monaco-editor/esm/vs/editor/contrib/find/findController.js';
// import 'monaco-editor/esm/vs/editor/contrib/folding/folding.js';
// import 'monaco-editor/esm/vs/editor/contrib/hover/hover.js';
// import 'monaco-editor/esm/vs/editor/contrib/multicursor/multicursor.js';
// import 'monaco-editor/esm/vs/editor/contrib/smartSelect/smartSelect.js';
// import 'monaco-editor/esm/vs/editor/contrib/snippet/snippetController2.js';
// import 'monaco-editor/esm/vs/editor/contrib/suggest/suggestController.js';
// import 'monaco-editor/esm/vs/editor/standalone/browser/quickOpen/quickCommand.js';
import * as monaco from 'monaco-editor/esm/vs/editor/editor.api.js';

// --------—--------—--------—--------—--------—--------—--------—
// CodeMirror
// --------—--------—--------—--------—--------—--------—--------—
// import { EditorState, EditorView, basicSetup } from "@codemirror/basic-setup"
// import { python } from "@codemirror/lang-python"
// onMount(() => {
// 	let editor = new EditorView({
// 	state: EditorState.create({
// 		extensions: [basicSetup, python()]
// 	}),
// 	parent: document.body
// 	})
// });

let containerElt;
onMount(() => {
	const editor = monaco.editor.create(containerElt, {
		value: `def test():\n\tprint(123)\n`,
		theme: "vs-light",
		language: "python",
		minimap: { enabled: false }
    });

	editor.addAction({
		id: 'my-unique-id',
		label: 'My Label!!!',
		keybindings: [
			monaco.KeyMod.CtrlCmd | monaco.KeyCode.Enter,
		],
		// A precondition for this action.
		precondition: null,
		// A rule to evaluate on top of the precondition in order to dispatch the keybindings.
		keybindingContext: null,
		contextMenuGroupId: 'navigation',
		contextMenuOrder: 1.5,
		// Method that will be executed when the action is triggered.
		// @param editor The editor instance is passed in as a convinience
		run: function(ed) {
			alert("i'm running => " + ed.getPosition());
			return null;
		}
	});


  });

</script>

<div bind:this={containerElt} id="container-editor" style="height: 80vh"></div>
