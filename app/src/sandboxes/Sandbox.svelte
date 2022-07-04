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
let flags = `-F \\t -v abc=2`; // ["-F", "\\t", "-v", "abc=2"];
let error;

// Supported flags
const FLAGS = [
	{
		name: "Set delimiter",
		flag: "-F",  // 
		options: [
			{ name: "Tabs", value: "\\t" },
			{ name: "Commas", value: "," },
			{ name: "Spaces", value: `" "` }
		]
	},
	{
		name: "Define Variable",
		flag: "-v",  // varname=value
		options: [
			{ name: "Add new variable", value: "myvar=123" },
		],
		multiple: true
	},
];

// Reactive logic
$: langCmd = tool === "jq" ? "json" : "cpp";
$: langIO = tool === "jq" ? "json" : null;
$: if(CLI.ready && input && command && tool && flags && $sandbox.settings.interactive) run();


// =============================================================================
// 
// =============================================================================

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
	if(busy)
		return;
	busy = true;

	// Prepare inputs
	const toolName = tool === "awk" ? "gawk" : tool;
	const params = [];
	if(tool === "jq")
		params.push("-M");
	params.push(command.trim());

const flagsArr = parseFlags(flags);
console.log("FLAGS", flagsArr);

	// Run 
	try {
		await CLI.fs.writeFile("sandbox", input);
console.warn([...flagsArr, ...params, "sandbox"]);
		const { stdout, stderr } = await CLI.exec(toolName, [...flagsArr, ...params, "sandbox"]);
		output = stdout;
		error = stderr;
	} catch (error) {
		console.error(error);
	} finally {
		busy = false;
	}
}

// Convert flags strings into array. Don't just do `.split(" ")` because a flag could
// be `-F " "`, which we want to treat as [`-F`, `" "`], not [`-F`, `"`, `"`]
function parseFlags(flags) {
	return flags.match(/[A-Za-z0-9-_=\\,\.]+|"[^"]+"/g);
	// // AST parser doesn't support equal sign in bash yet
	// str = str.replaceAll("=", "{EQUAL}");

	// // Parse as AST and retrieve array of arguments
	// const args = parse(`awk ${str} --end`)[0].args.map(d => d.value.replaceAll("{EQUAL}", "="));
	// console.warn("args", args);
	// return args;
}

//
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
						{#each FLAGS as flag, i}
							<DropdownItem header>{flag.name}</DropdownItem>
							{#each flag.options || [] as option}
								<DropdownItem on:click={() => {
									setFlag(flag, option.value)
								}}>{option.name}</DropdownItem>
							{/each}
							{#if i < FLAGS.length - 1}
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
