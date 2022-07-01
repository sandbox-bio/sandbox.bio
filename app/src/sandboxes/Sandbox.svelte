<script>
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";
import { Button, ButtonDropdown, DropdownItem, DropdownMenu, DropdownToggle, Input, Tooltip } from "sveltestrap";
import { sandbox } from "stores/sandbox";
import IDE from "components/IDE.svelte";
import awk_data from "./orders.txt";

export let tool = "jq";

// Tools to load in playground
const TOOLS = [
	{ loading: "eager", tool: "jq", version: "1.6" },
	{ loading: "eager", tool: "gawk", version: "5.1.0", reinit: true },
	{ loading: "eager", tool: "grep", version: "3.7", reinit: true  },
	{ loading: "eager", tool: "sed", version: "4.8", reinit: true  },
];

// State
let CLI = {};
let busy = false;
let divSettingEnter;

let command = `/Burrito/ { print $3 }`;
let input = awk_data;
let output;
let flags = ["-F", "\\t", "-v", "abc=2"];
let error;


const FLAGS  =[
	{
		name: "Set delimiter",
		flags: "-F",  // 
		values: [
			{ name: "Spaces", value: " " },
			{ name: "Tabs", value: "\t" },
			{ name: "Commas", value: "," }
		]
	},
	{
		name: "Add Variable",
		flag: "-v",  // varname=value
		multiple: true
	},
];

// Reactive logic
$: langCmd = tool === "jq" ? "json" : "cpp";
$: langIO = tool === "jq" ? "json" : null;
$: if(CLI.ready && input && command && tool && flags && $sandbox.settings.interactive) run();


// Initialize sandbox
onMount(async () => {
	// Initialize store with data from localforage
	await sandbox.init();

	// Initialize Aioli
	CLI = await new Aioli(TOOLS, {
		env: ["localhost", "dev.sandbox.bio"].includes(window.location.hostname) ? "stg" : "prd",
		printInterleaved: false
	});
	CLI.ready = true;
});

// Run command with given input and show resulting output/error
async function run() {
	busy = true;

	// Prepare inputs
	const toolName = tool === "awk" ? "gawk" : tool;
	const params = [];
	if(tool === "jq")
		params.push("-M");
	params.push(command.trim());

	// Run 
	try {
		await CLI.fs.writeFile("sandbox", input);
		const { stdout, stderr } = await CLI.exec(toolName, [...flags, ...params, "sandbox"]);
		output = stdout;
		error = stderr;
	} catch (error) {
		console.error(error);
	} finally {
		busy = false;
	}
}
</script>

<h4>{tool} sandbox</h4>

<div class="row">
	<!-- Command -->
	<div class="col-md-6">
		<div class="row ide mb-4 mt-4">
			<div class="d-flex flex-row">
				<div class="pe-4">
					<h5>Command</h5>
				</div>
				<div bind:this={divSettingEnter} class="pe-4">
					<Input type="switch" label="Interactive" bind:checked={$sandbox.settings.interactive} />
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
					<Button color="primary" size="sm" on:click={run} disabled={busy}>
						Run
					</Button>
				</div>
			</div>

			<!-- Errors -->
			{#if error}
				<pre class="text-danger pre-scrollable">{error}</pre>
			{/if}
		</div>
	</div>

	<!-- Flags -->
	<div class="col-md-6">
		<div class="row ide mb-4 mt-4">
			<div class="d-flex flex-row">
				<div class="pe-4">
					<h5>Flags</h5>
				</div>
				<ButtonDropdown size="sm">
					<DropdownToggle color="primary" caret>Add flag</DropdownToggle>
					<DropdownMenu>
						{#each FLAGS as flag}
							<DropdownItem>{flag.name}</DropdownItem>
						{/each}
					</DropdownMenu>
				</ButtonDropdown>
			</div>

			<IDE
				lang={null}
				code={flags.join(" ")}
				on:update={d => flags = d.detail.split(" ")} />
		</div>
	</div>
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
</style>
