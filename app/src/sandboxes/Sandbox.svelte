<script>
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";
import { Button, ButtonDropdown, DropdownItem, DropdownMenu, DropdownToggle, Input, Tooltip } from "sveltestrap";
import { tool, data, sandbox, TOOLS, FLAGS, FLAG_SETTING, FLAG_BOOLEAN, FLAG_PARAM } from "stores/sandbox";
import IDE from "components/IDE.svelte";

// State
let CLI = {};
let busy = true;
let divSettingEnter;
let output;
let error;

// Reactive logic
$: langCmd = $tool?.name === "jq" ? "json" : "cpp";
$: langIO = $tool?.name === "jq" ? "json" : null;
$: if(CLI.ready && $data.input && $data.command && $tool?.name && $sandbox.settings.interactive) run($data.flags);


// =============================================================================
// Main logic
// =============================================================================

// Initialize sandbox
onMount(async () => {
	// Initialize store with data from localforage
	await sandbox.init();

	// Get tool from URL
	const params = new URL(window.location).searchParams;
	$tool = TOOLS.find(d => d.name === (params.get("id") || "awk"));
	if(!$tool)
		throw "Unrecognized tool";
	$sandbox = $sandbox;  // force update the $data.flags derived store

	// Initialize Aioli
	CLI = await new Aioli($tool.aioliConfig, {
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
	if($tool.name === "jq")
		params.push("-M");
	// Add user flags
	params = params.concat(parseFlags($data.flags)).map(d => d.replaceAll('"', ''));
	// Add user command
	params.push($data.command.trim());
	// Add file to operate on
	params.push("sandbox");

	// Run
	try {
		await CLI.fs.writeFile("sandbox", $data.input);
		const { stdout, stderr } = await CLI.exec($tool.aioliConfig.tool, params);
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

function getFlag(option) {
	const flagsArr = parseFlags($data.flags);
	const flagIndex = flagsArr.findIndex(d => d === option.flag);

	// Boolean flags: only true if the flag is present
	if(option.type === FLAG_BOOLEAN)
		return flagIndex !== -1;

	// Setting flags: return the value if the flag is present
	if(option.type === FLAG_SETTING)
		return flagIndex !== -1 ? flagsArr[flagIndex + 1] : null;

	throw "Should only use getFlag with option types 'boolean' or 'setting'";
}

// Add or modify a flag (`value` only used by "setting" flag)
function setFlag(option, value) {
	let flagsArr = parseFlags($data.flags);
	const flagIndex = flagsArr.findIndex(d => d === option.flag);

	// Boolean flags: toggle true/false
	if(option.type === FLAG_BOOLEAN) {
		if(flagIndex !== -1)
			delete flagsArr[flagIndex];
		else
			flagsArr.push(option.flag);

	// Setting flag: add or change existing value
	} else if(option.type === FLAG_SETTING) {
		if(flagIndex !== -1)
			flagsArr[flagIndex + 1] = value;
		else
			flagsArr = flagsArr.concat([ option.flag, value ]);

	// Param flag: add a new param
	} else if(option.type === FLAG_PARAM) {
		flagsArr = flagsArr.concat([ option.flag, option.value ]);
	}

	$sandbox.data[$tool.name].flags = flagsArr.join(" ").trim();
}


// =============================================================================
// HTML
// =============================================================================
</script>

{#if $tool}
	<h4 class="mb-0">{$tool.name} sandbox</h4>

	<div class="row">
		<!-- Command -->
		<div class="col-md-6">
			<div class="row ide mb-4 mt-4">
				<div class="d-flex flex-row mb-2">
					<div class="pe-3 pt-1 pb-1">
						<h5>Command</h5>
					</div>
					<div bind:this={divSettingEnter} class="pt-1">
						<Input type="checkbox" label="Interactive" bind:checked={$sandbox.settings.interactive} />
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
							code={$data.command}
							on:update={d => $sandbox.data[$tool.name].command = d.detail}
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
					<div class="pe-1 pt-2">
						<h5>Flags</h5>
					</div>
					{#each FLAGS[$tool.name] as option}
						<!-- Setting -->
						{#if option.type === FLAG_SETTING}
							<ButtonDropdown size="sm" class="mx-1 my-1">
								<DropdownToggle color="primary" caret>{option.name}</DropdownToggle>
								<DropdownMenu>
									{#each option.values as value}
										<DropdownItem on:click={() => setFlag(option, value.value)} class="small">
											{value.name}: <code>{value.value}</code>
										</DropdownItem>
									{/each}
								</DropdownMenu>
							</ButtonDropdown>

						<!-- Boolean Setting -->
						{:else if option.type === FLAG_BOOLEAN}
							<div class="mx-1 pt-2">
								<Input type="checkbox" label={option.name} checked={getFlag(option)} on:change={() => setFlag(option)} />
							</div>	

						<!-- Parameters -->
						{:else if option.type === FLAG_PARAM}
							<Button size="sm" class="mx-1 my-1" on:click={() => setFlag(option)}>+ {option.name}</Button>
						{/if}
					{/each}
				</div>

				<IDE
					lang={null}
					code={$data.flags}
					on:update={d => $sandbox.data[$tool.name].flags = d.detail} />
			</div>
		</div>
	</div>

	<!-- Input / Output -->
	<div class="row">
		<div class="col-md-6 ide">
			<h5>Input</h5>
			<IDE
				lang={langIO}
				code={$data.input}
				on:update={d => $sandbox.data[$tool.name].input = d.detail} />
		</div>
		<div class="col-md-6 ide">
			<h5>Output</h5>
			<IDE
				lang={langIO}
				code={output}
				on:update={d => output = d.detail} editable={false} />
		</div>
	</div>
{/if}

<style>
.ide {
	font-size: 15px;  /* default = 16px */
}
</style>
