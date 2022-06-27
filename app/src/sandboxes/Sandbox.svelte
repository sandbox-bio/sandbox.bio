<script>
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";
import { Button, Input, Tooltip } from "sveltestrap";
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
let settings = {
	interactive: true
};

let command = `. | length\n`;
let input = JSON.stringify({"abc": "sdf", "def": "fsd", "ghi": { "a": 123, "b": 456 }}, null, 2);
let output;
let error;

let divSettingEnter;

// Reactive logic
$: langCmd = tool === "jq" ? "json" : "cpp";
$: langIO = tool === "jq" ? "json" : null;
$: if(CLI.ready && input && command && tool && settings.interactive) run();

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
	const params = [];
	if(tool === "jq")
		params.push("-M");

	// Prepare inputs
	const commandTxt = command.trim();

	// Run 
	try {
		await CLI.fs.writeFile("sandbox", input);
		const { stdout, stderr } = await CLI.exec(tool, [...params, commandTxt, "sandbox"]);
		output = stdout;
		error = stderr;
	} catch (error) {
		console.error(error);
	}
}
</script>

<!-- Select a tool -->
<div class="form-floating col-md-2">
	<select class="form-select" id="tool" bind:value={tool}>
		<option value="jq">jq</option>
		<option value="gawk">awk</option>
		<option value="grep">grep</option>
		<option value="sed">sed</option>
	</select>
	<label for="tool">Choose a tool</label>
</div>

<!-- Command -->
<div class="row ide ide-command mb-4 mt-4">
	<div class="d-flex flex-row">
		<div class="pe-4">
			<h5>Command</h5>
		</div>
		<div bind:this={divSettingEnter}>
			<Input type="switch" label="Interactive" bind:checked={settings.interactive} />
		</div>
		<Tooltip target={divSettingEnter}>
			Run after each keypress
		</Tooltip>
	</div>

	<!-- Command box -->
	<div class="d-flex flex-row">
		<div class="w-100">
			<IDE
				lang={langCmd}
				code={command}
				on:update={d => command = d.detail}
				on:run={run} />
		</div>
		<div class="flex-shrink-1 ps-3">
			<Button color="primary" size="lg" on:click={run}>Run</Button>
		</div>
	</div>

	<!-- Errors -->
	{#if error}
		<pre class="text-danger pre-scrollable">{error}</pre>
	{/if}
</div>

<!-- Input / Output -->
<div class="row">
	<div class="col-md-6 ide">
		<h5>Input</h5>
		<IDE
			lang={langIO}
			code={input}
			on:update={d => input = d.detail} />
	</div>
	<div class="col-md-6 ide">
		<h5>Output</h5>
		<IDE
			lang={langIO}
			code={output}
			on:update={d => output = d.detail} editable={false} />
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
