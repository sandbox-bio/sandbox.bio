<script>
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";
import { Button, ButtonDropdown, DropdownItem, DropdownMenu, DropdownToggle, Input, Spinner, Tooltip } from "sveltestrap";
import { tool, data, sandbox, TOOLS, EXAMPLES, FLAGS, FLAG_SETTING, FLAG_BOOLEAN, FLAG_PARAM } from "$stores/sandbox";
import { config } from "$stores/config";
import IDE from "$components/IDE.svelte";

// State
let CLI = {};
let ready = false;    // whether CLI is done loading
let busy = false;     // whether CLI is busy
let busy_ui = false;  // whether UI is busy (CLI.busy is whether Aioli is busy)
let divSettingEnter;
let output;
let error;

// Reactive logic
$: langCmd = { awk: "awk", jq: "json" }[$tool?.name];
$: langIO = $tool?.name === "jq" ? "json" : null;
$: if(ready && $data.input && $data.command !== null && $tool?.name && $sandbox.settings.interactive) run($data.flags);

// If update flags from input box, need to update checkboxes!
$: if($data?.flags !== null) {
	$tool = $tool;
}

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
	CLI = await new Aioli([
		{ tool: "base", version: "1.0.0" },	
		$tool.aioliConfig
	], {
		env: ["localhost", "dev.sandbox.bio"].includes(window.location.hostname) ? "stg" : "prd",
		printInterleaved: false
	});
	ready = true;
	busy = false;
});

// Run command with given input and show resulting output/error
async function run() {
	if(!ready || busy)
		return;

	// If the CLI is still busy after some time, show a spinner, otherwise don't (avoids flickering)
	busy = true;
	setTimeout(() => {
		if(busy)
			busy_ui = true;
	}, 300);

	// Prepare parameters
	let params = [];
	if($tool.name === "jq")
		params.push("-M");
	// Add user flags
	params = params.concat(parseFlags($data.flags)).map(d => d.replaceAll('"', ''));
	// Add user command
	if($data.command.trim())
		params.push($data.command.trim());
	// Add file to operate on
	params.push("sandbox");

	// Run
	try {
		await CLI.fs.writeFile("sandbox", $data.input);
		const { stdout, stderr } = await CLI.exec($tool.aioliConfig.tool, params);
		output = stdout;
		error = stderr;

		// Analytics
		try {
			const d = {
				playground: $tool.name,
				example: EXAMPLES[$tool.name].some(example => example.input == $data.input && example.flags == $data.flags && example.command == $data.command)
			};
			fetch(`${$config.api}/ping`, {
				method: "POST",
				mode: "no-cors",
				body: JSON.stringify(d)
			});
		} catch (error) {}
	} catch (error) {
		console.error(error);
	} finally {
		busy = false;
		busy_ui = false;
	}
}

function updateVar(varName, value) {
	$sandbox.data[$tool.name][varName] = value;
}


// =============================================================================
// Flag Management
// =============================================================================

// Convert flags strings into array. Don't just do `.split(" ")` because a flag could
// be `-F " "`, which we want to treat as [`-F`, `" "`], not [`-F`, `"`, `"`]
function parseFlags(flags) {
	// Source: https://stackoverflow.com/a/16261693
	// Note that the AST parser doesn't support equal sign in bash yet
	// and "shell-quote"'s ShellQuote.parse() doesn't process \t properly
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

	updateVar("flags", flagsArr.join(" ").trim());
}


// =============================================================================
// HTML
// =============================================================================
</script>

{#if $tool}
	<h4 class="mb-0">
		{$tool.name} sandbox

		<ButtonDropdown class="mx-1 my-1">
			<DropdownToggle color="primary" caret>Examples</DropdownToggle>
			<DropdownMenu>
				{#each EXAMPLES[$tool.name] as example}
					{@const active = example.input == $data.input && example.flags == $data.flags && example.command == $data.command}
					<DropdownItem class="py-2" active={active} on:click={() => {
						updateVar("input", example.input);
						updateVar("flags", example.flags);
						updateVar("command", example.command);
					}}>
						{example.name}
					</DropdownItem>
				{/each}
			</DropdownMenu>
		</ButtonDropdown>
	</h4>

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
				<div class="d-flex w-100" style="border:0px solid red">
					<div class="col-11" style="border:0px solid blue">
						<IDE
							lang={langCmd}
							code={$data.command}
							on:update={d => updateVar("command", d.detail)}
							on:run={run} />
					</div>
					<div class="col-1" style="border:0px solid green">
						{#if $sandbox.settings.interactive}
							{#if busy_ui || !ready}
								<Spinner class="ms-3" color="primary" />
							{/if}
						{:else}
							<Button color="primary" size="sm" on:click={run} disabled={busy_ui}>
								Run
							</Button>
						{/if}
					</div>
				</div>				  

				<!-- Errors -->
				{#if error}
					<pre class="text-danger">{error}</pre>
				{/if}
			</div>

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

			<div class="ide">
				<h5>Input</h5>
				<IDE
					lang={langIO}
					code={$data.input}
					on:update={d => updateVar("input", d.detail)} />
			</div>	
		</div>

		<!-- Flags -->
		<div class="col-md-6">
			<div class="ide mt-4">
				<h5>Output</h5>
				<IDE
					lang={langIO}
					code={output}
					on:update={d => output = d.detail} editable={false} />
			</div>
		</div>
	</div>
{/if}

<style>
.ide {
	font-size: 15px;  /* default = 16px */
}

pre {
	white-space: pre-wrap;  /* to avoid scrolling horizontally for long error messages */
}
</style>
