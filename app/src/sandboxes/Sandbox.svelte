<script>
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";
import { Button, ButtonDropdown, DropdownItem, DropdownMenu, DropdownToggle, Input, Tooltip } from "sveltestrap";
import { sandbox, TOOLS, FLAGS } from "stores/sandbox";
import IDE from "components/IDE.svelte";
import awk_data from "./orders.txt";

export let tool = "awk";

// State
let CLI = {};
let busy = true;
let divSettingEnter;

let command = `/Burrito/ { print abc, $3 }`;  // `.`;
let input = awk_data;  // `{"a": 4, "b":      5}`;
let flags = `-F \\t -v abc="I like"`;  // '';
let output;
let error;

// Reactive logic
$: toolAioli = TOOLS.find(d => d.name === tool)?.aioliConfig;
$: langCmd = tool === "jq" ? "json" : "cpp";
$: langIO = tool === "jq" ? "json" : null;
$: if(CLI.ready && input && command && tool && $sandbox.settings.interactive) run(flags);


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
	params = params.concat(parseFlags(flags)).map(d => d.replaceAll('"', ''));
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
	// Source: https://stackoverflow.com/a/16261693
	// Note that the AST parser doesn't support equal sign in bash yet
	return flags.match(/(?:[^\s"]+|"[^"]*")+/g) || [];
}

// Add or modify a flag
async function setFlag(option) {
	let flagsArr = parseFlags(flags);

	// Only allow a single copy of this flag
	if(!option.multiple) {
		const flagIndex = flagsArr.findIndex(d => d === option.flag);
		// If flag doesn't exist yet
		if(flagIndex === -1) {
			if(option.boolean) flagsArr.push(option.flag);
			else flagsArr = flagsArr.concat([ option.flag, option.value ]);
		// If flag exists, overwrite its value
		} else {
			if(option.boolean) delete flagsArr[flagIndex];
			else flagsArr[flagIndex + 1] = option.value;
		}

	// Allow multiple copies of this flag (e.g. `-v` in awk for defining variables)
	} else {
		flagsArr = flagsArr.concat([ option.flag, option.value ]);
	}

	flags = flagsArr.join(" ").trim();
}


// =============================================================================
// HTML
// =============================================================================
</script>

<h4 class="mb-0">{tool} sandbox</h4>

<div class="row">
	<!-- Command -->
	<div class="col-md-6">
		<div class="row ide mb-4 mt-4">
			<div class="d-flex flex-row mb-2">
				<div class="pe-3 pt-1 pb-1">
					<h5>Command</h5>
				</div>
				<div bind:this={divSettingEnter} class="pt-1">
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
			<div class="d-flex flex-row mb-2">
				<div class="pe-2 pt-2">
					<h5>Flags</h5>
				</div>
				<ButtonDropdown size="sm">
					<DropdownToggle color="primary" caret>{tool} flags</DropdownToggle>
					<DropdownMenu>
						{#each FLAGS[tool] as category, i}
							<DropdownItem header class="text-primary">{category.name}</DropdownItem>
							{#each category.options || [] as option}
								<DropdownItem on:click={() => setFlag(option)}>
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
