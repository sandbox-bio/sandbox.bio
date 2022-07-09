<script>
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";
import { Button, ButtonDropdown, DropdownItem, DropdownMenu, DropdownToggle, Input, Tooltip } from "sveltestrap";
import { sandbox } from "stores/sandbox";
import IDE from "components/IDE.svelte";
import awk_data from "./orders.txt";

export let tool = "awk";

// Tools to load in playground
const TOOLS = [
	{ name: "jq", aioliConfig: { tool: "jq", version: "1.6" }},
	{ name: "awk", aioliConfig: { tool: "gawk", version: "5.1.0", reinit: true }},
	{ name: "grep", aioliConfig: { tool: "grep", version: "3.7", reinit: true }},
	{ name: "sed", aioliConfig: { tool: "sed", version: "4.8", reinit: true }}
];

// State
let CLI = {};
let busy = true;
let divSettingEnter;

let command = `/Burrito/ { print $3 }`;
let input = awk_data;
let output;
let flags = `-F \\t -v abc=2`;
let error;

// Supported flags
const FLAGS = {
	awk: [
		{
			name: "Set delimiter",
			flag: "-F",
			options: [
				{ name: "Tabs", value: "\\t" },
				{ name: "Commas", value: "," },
				// { name: "Spaces", value: `" "` }
			]
		},
		{
			name: "Define Variable",
			flag: "-v",
			options: [
				{ name: "Add new variable", value: "myvar=123" },
			],
			multiple: true
		}
	]
}

// Reactive logic
$: toolAioli = TOOLS.find(d => d.name === tool)?.aioliConfig;
$: langCmd = tool === "jq" ? "json" : "cpp";
$: langIO = tool === "jq" ? "json" : null;
$: if(CLI.ready && input && command && tool && flags && $sandbox.settings.interactive) run();


// =============================================================================
// Main logic
// =============================================================================

// Initialize sandbox
onMount(async () => {
	// Initialize store with data from localforage
	await sandbox.init();

	// Initialize Aioli
	CLI = await new Aioli(toolAioli, {
		env: ["localhost", "dev.sandbox.bio"].includes(window.location.hostname) ? "stg" : "prd",
		printInterleaved: false
	});
	CLI.ready = true;
	busy = false;
});

// Run command with given input and show resulting output/error
async function run() {
	if(busy)
		return;
	busy = true;

	// Prepare parameters
	let params = [];
	if(tool === "jq")
		params.push("-M");
	// Add user flags
	params = params.concat(parseFlags(flags));
	// Add user command
	params.push(command.trim());
	// Add file to operate on
	params.push("sandbox");

	// Run
	try {
		await CLI.fs.writeFile("sandbox", input);
		const { stdout, stderr } = await CLI.exec(toolAioli.tool, params);
		output = stdout;
		error = stderr;
	} catch (error) {
		console.error(error);
	} finally {
		busy = false;
	}
}


// =============================================================================
// Flag Management
// =============================================================================

// Convert flags strings into array. Don't just do `.split(" ")` because a flag could
// be `-F " "`, which we want to treat as [`-F`, `" "`], not [`-F`, `"`, `"`]
function parseFlags(flags) {
	// Note that the AST parser doesn't support equal sign in bash yet
	return flags.match(/[A-Za-z0-9-_=\\,\.]+|"[^"]+"/g);
}

// Add or modify a flag
async function setFlag(flag, value) {
	console.log("flag", flag, value)
	let flagsArr = parseFlags(flags);

	// Only allow a single copy of this flag
	if(!flag.multiple) {
		const flagIndex = flagsArr.findIndex(d => d === flag.flag);
		if(flagIndex === -1)
			return;
		console.log(flagIndex)

		if(flag.boolean) {
			// TODO:
		} else {
			flagsArr[flagIndex + 1] = value;
		}
		console.log(flagsArr);

	// Allow multiple copies of this flag (e.g. `-v` in awk for defining variables)
	} else {
		flagsArr.push(flag.flag);
		flagsArr.push(value);
	}

	console.warn("flagsArr", flagsArr);
	flags = flagsArr.join(" ");
}


// =============================================================================
// HTML
// =============================================================================
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
					<DropdownToggle color="primary" caret>{tool} flags</DropdownToggle>
					<DropdownMenu>
						{#each FLAGS[tool] as flag, i}
							<DropdownItem header class="text-primary">{flag.name}</DropdownItem>
							{#each flag.options || [] as option}
								<DropdownItem on:click={() => setFlag(flag, option.value)}>
									&bullet; {option.name}
								</DropdownItem>
							{/each}
							{#if i < FLAGS[tool].length - 1}
								<DropdownItem divider />
							{/if}
						{/each}
					</DropdownMenu>
				</ButtonDropdown>
			</div>

			<IDE
				lang={null}
				code={flags}
				on:update={d => flags = d.detail} />
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
