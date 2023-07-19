<script>
import { Button, ButtonDropdown, DropdownItem, DropdownMenu, DropdownToggle, Form, FormGroup, Input, Spinner } from "sveltestrap";
import IDE from "$components/IDE.svelte";
import Setting from "$components/playgrounds/Setting.svelte";
import { sandbox, EXAMPLES, FLAGS, FLAG_SETTING, FLAG_BOOLEAN, FLAG_PARAM } from "$stores/sandbox";

export let tool = {}; // tool info, defined in stores/sandbox.js
export let CLI = {}; // initialized Aioli object
export let ready = false; // whether Aioli is ready to go

// State
let busy = false; // whether CLI is busy
let busy_ui = false; // whether UI is busy (CLI.busy is whether Aioli is busy)
let output;
let error;

// Reactive logic
$: langCmd = { awk: "awk", jq: "json" }[tool.name];
$: langIO = tool.name === "jq" ? "json" : null;
$: userInput = $sandbox.data[tool.name];
$: flags = parseFlags(userInput.flags);
$: hasExample = EXAMPLES[tool.name].findIndex(
	(example) => example.input === userInput.input && example.flags === userInput.flags && example.command === userInput.command
);

// Re-run if any user input changes. Need the "if ready" to run command when Aioli
// is loaded, otherwise it won't run anything on first load.
$: if (ready && $sandbox.interactive && userInput) run();

// =============================================================================
// Main logic
// =============================================================================

// Run command with given input and show resulting output/error
async function run() {
	// Don't run a command if we're already running another one
	if (busy) return;

	// If the CLI is still busy after some time, show a spinner, otherwise don't (avoids flicker)
	busy = true;
	setTimeout(() => {
		if (busy) busy_ui = true;
	}, 300);

	// Prepare parameters
	let params = [];
	if (tool.name === "jq") params.push("-M");
	// Add user flags
	params = params.concat(flags).map((d) => d.replaceAll('"', ""));
	// Add user command
	if (userInput.command.trim()) params.push(userInput.command.trim());
	// Add file to operate on
	params.push("sandbox");

	// Run
	try {
		await CLI.fs.writeFile("sandbox", userInput.input);
		const { stdout, stderr } = await CLI.exec(tool.aioli.tool, params);
		output = stdout;
		error = stderr;

		// Analytics
		try {
			const d = {
				playground: tool.name,
				example: hasExample > -1
			};
			fetch(`/ping`, {
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

// Updating `userInput.command` doesn't save the changes to localforage, so we
// need to explicitely save them.
function updateUserInput(varName, value) {
	$sandbox.data[tool.name][varName] = value;
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
	const flagIndex = flags.findIndex((d) => d === option.flag);

	// Boolean flags: only true if the flag is present
	if (option.type === FLAG_BOOLEAN) return flagIndex !== -1;

	// Setting flags: return the value if the flag is present
	if (option.type === FLAG_SETTING) return flagIndex !== -1 ? flags[flagIndex + 1] : null;

	throw "Should only use getFlag with option types 'boolean' or 'setting'";
}

// Add or modify a flag (`value` only used by "setting" flag)
function setFlag(option, value) {
	const flagIndex = flags.findIndex((d) => d === option.flag);

	// Boolean flags: toggle true/false
	if (option.type === FLAG_BOOLEAN) {
		if (flagIndex !== -1) delete flags[flagIndex];
		else flags.push(option.flag);

		// Setting flag: add or change existing value
	} else if (option.type === FLAG_SETTING) {
		if (flagIndex !== -1) flags[flagIndex + 1] = value;
		else flags = flags.concat([option.flag, value]);

		// Param flag: add a new param
	} else if (option.type === FLAG_PARAM) {
		flags = flags.concat([option.flag, option.value]);
	}

	updateUserInput("flags", flags.join(" ").trim());
}

// =============================================================================
// HTML
// =============================================================================
</script>

<!-- Title and links to other playgrounds -->
<div class="d-flex">
	<div class="d-flex">
		<h4 class="mb-0">{tool.name} sandbox</h4>
	</div>
	<div class="d-flex ms-auto">
		<slot name="playgrounds" />
	</div>
</div>

<!-- Examples -->
<div class="row g-2 mt-1">
	<div class="col-auto">
		<label for="examples" class="col-form-label text-muted">Examples:</label>
	</div>
	<div class="col-auto">
		<Input
			id="examples"
			type="select"
			size="sm"
			value={hasExample}
			on:change={(d) => {
				const example = EXAMPLES[tool.name][+d.target.value];
				updateUserInput("input", example.input);
				updateUserInput("flags", example.flags);
				updateUserInput("command", example.command);
			}}
		>
			{#each EXAMPLES[tool.name] as example, i}
				<option value={i}>{example.name}</option>
			{/each}
		</Input>
	</div>
</div>

<!-- Playground -->
<div class="row">
	<!-- Input -->
	<div class="col-md-6">
		<!-- Command -->
		<div class="row ide mb-4 mt-4">
			<div class="d-flex flex-row mb-2">
				<div class="pe-3 pt-1 pb-1">
					<h5>Command</h5>
				</div>
				<Setting tooltip="Run after each keypress">
					<Input type="checkbox" label="Interactive" bind:checked={$sandbox.interactive} />
				</Setting>
			</div>

			<!-- Command box -->
			<div class="d-flex w-100">
				<div class="col-11">
					<IDE lang={langCmd} code={userInput.command} on:update={(d) => updateUserInput("command", d.detail)} on:run={run} />
				</div>
				<div class="col-1">
					{#if $sandbox.interactive}
						{#if busy_ui}
							<Spinner class="ms-3" color="primary" />
						{/if}
					{:else}
						<Button color="primary" size="sm" on:click={run} disabled={busy_ui}>Run</Button>
					{/if}
				</div>
			</div>

			<!-- Errors -->
			{#if error}
				<pre class="text-danger">{error}</pre>
			{/if}
		</div>

		<!-- Flags -->
		<div class="row ide mb-4 mt-4">
			<div class="d-flex flex-row mb-2">
				<div class="pe-1 pt-2">
					<h5>Flags</h5>
				</div>
				{#key userInput.flags}
					{#each FLAGS[tool.name] as option}
						<!-- Setting -->
						{#if option.type === FLAG_SETTING}
							<ButtonDropdown size="sm" class="mx-1 my-1">
								<DropdownToggle color="outline-primary" caret>{option.name}</DropdownToggle>
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
							<div class="mx-2 pt-1">
								<Setting tooltip={option.description}>
									<Input type="checkbox" label={option.name} checked={getFlag(option)} on:change={() => setFlag(option)} />
								</Setting>
							</div>

							<!-- Parameters -->
						{:else if option.type === FLAG_PARAM}
							<Button size="sm" class="mx-1 my-1" on:click={() => setFlag(option)}>+ {option.name}</Button>
						{/if}
					{/each}
				{/key}
			</div>

			<IDE lang={null} code={userInput.flags} on:update={(d) => ($sandbox.data[tool.name].flags = d.detail)} />
		</div>

		<div class="ide">
			<h5>Input</h5>
			<IDE lang={langIO} code={userInput.input} on:update={(d) => updateUserInput("input", d.detail)} />
		</div>
	</div>

	<!-- Result -->
	<div class="col-md-6">
		<div class="ide mt-4">
			<h5>Output</h5>
			<IDE lang={langIO} code={output} on:update={(d) => (output = d.detail)} editable={false} />
		</div>
	</div>
</div>

<style>
.ide {
	font-size: 15px; /* default = 16px */
}

pre {
	white-space: pre-wrap; /* to avoid scrolling horizontally for long error messages */
}
</style>
