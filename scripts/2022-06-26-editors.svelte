<script>
import { onMount } from 'svelte';

// // =============================================================================
// // Highlight.js
// import hljs from 'highlight.js/lib/core';
// import json from 'highlight.js/lib/languages/json';
// import awk from 'highlight.js/lib/languages/awk';
// import 'highlight.js/styles/github.css';
// hljs.registerLanguage('json', json);
// hljs.registerLanguage('awk', awk);
// onMount(() => {
// 	hljs.highlightAll();
// })


// =============================================================================
// CodeMirror
// npm dedupe
// Can't find good documentation for CodeMirror... + Won't initialize
//    Uncaught Error: Unrecognized extension value in extension set ([object Object]). This sometimes happens because multiple instances of @codemirror/state are loaded, breaking instanceof checks
// import { EditorState, EditorView, basicSetup } from "@codemirror/basic-setup"
// import {syntaxHighlighting, defaultHighlightStyle} from "@codemirror/language"
import { EditorView, basicSetup, minimalSetup } from 'codemirror'
import { EditorState } from '@codemirror/state'
// import { lineNumbers } from '@codemirror/view'

import { json, jsonParseLinter } from "@codemirror/lang-json";
let elCodeMirror;
onMount(() => {

	// https://github.com/codemirror/basic-setup/blob/main/src/codemirror.ts#L47

	let editor = new EditorView({
		state: EditorState.create({
			extensions: [basicSetup.slice(2), json()],
			doc: `{
  "abc": "sdf",
  "bcd": "sdf",
  "csdf": {
    "a": 4,
    "b": "3"
  }
}`,
		}),
		parent: elCodeMirror //
		// extensions: minimalSetup,
		// extensions: [basicSetup, python()],
	});
	console.log(editor);
});




// =============================================================================
// Ace

</script>

<h3>Highlight.js</h3>
<pre><code class="language-json">
{"{"}
  "abc": [],
  "bd": 23
{"}"}
</code></pre>
<pre><code class="language-awk">
gsub("abc","def", $1)
</code></pre>
	
<h3>CodeMirror</h3>
<div bind:this={elCodeMirror} style="border:1px solid #ccc" />

<h3>Ace</h3>
