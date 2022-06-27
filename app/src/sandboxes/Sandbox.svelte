<script>
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";
import IDE from "../components/IDE.svelte";

// Tools to load in playground
const TOOLS = [
	{ loading: "eager", tool: "jq", version: "1.6" },
	{ loading: "eager", tool: "gawk", version: "5.1.0", reinit: true },
	{ loading: "eager", tool: "grep", version: "3.7", reinit: true  },
	{ loading: "eager", tool: "sed", version: "4.8", reinit: true  },
];

// State
let CLI = {};
let tool = "jq";

let command = `. | length`;
let input = JSON.stringify({"abc": "sdf", "def": "fsd", "ghi": { "a": 123, "b": 456 }}, null, 2);
let output;
let error;

// Reactive logic
$: langCmd = tool === "jq" ? "json" : "cpp";
$: langIO = tool === "jq" ? "json" : null;
$: if(CLI.ready && input && command && tool) run();

// Initialize Aioli
onMount(async () => {
	CLI = await new Aioli(TOOLS, {
		env: ["localhost", "dev.sandbox.bio"].includes(window.location.hostname) ? "stg" : "prd",
		printInterleaved: false
	});
	CLI.ready = true;
});

// Run command with given input and show resulting output/error
async function run() {
	await CLI.fs.writeFile("sandbox", input);

	const { stdout, stderr } = await CLI.exec(tool, ["-M", command, "sandbox"]);
	if(!stderr)
		output = stdout;
	error = stderr;
}
</script>

<div class="form-floating col-md-3">
	<select class="form-select" id="tool" bind:value={tool}>
		<option value="jq">jq</option>
		<option value="awk">awk</option>
		<option value="grep">grep</option>
		<option value="sed">sed</option>
	</select>
	<label for="tool">Choose a tool</label>
</div>

<div class="row ide ide-command mb-4 mt-4">
	<h5>Command</h5>
	<IDE lang={langCmd} code={command} on:update={d => command = d.detail} />
	{#if error}
		<p class="text-danger">{error}</p>
	{/if}
</div>

<div class="row">
	<div class="col-md-6 ide">
		<h5>Input</h5>
		<IDE lang={langIO} code={input} on:update={d => input = d.detail} />
	</div>
	<div class="col-md-6 ide">
		<h5>Output</h5>
		<IDE lang={langIO} code={output} on:update={d => output = d.detail} editable={false} />
	</div>
</div>

<style>
.ide {
	font-size: 15px;  /* default = 16px */
}

.ide-command {
	max-height: 120px;
	overflow: scroll;
}
</style>
