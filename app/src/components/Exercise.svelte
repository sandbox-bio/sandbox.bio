<script>
import { status } from "./stores/status";
import { CLI } from "terminal/cli";
import { FormGroup, Icon, Input, Label, ListGroup, ListGroupItem, Spinner } from "sveltestrap";

export let criteria = [];        // List of criteria that must be true for the exercise to be complete
export let hints = [];           // Hints to show user
let statuses = [];               // Status of each criteria (true/false)
let busy = false;                // Whether the user manually asked to check their work
let nbHints = 0;

// Check status every time a command finishes running
$: if($status.terminal == "execDone"){
	console.log("Checking solutions...")
	check();
};

$: isDone = statuses.filter(d => d).length == statuses.length && statuses.length != 0;

// Validate user's input
async function check(manual=false)
{
	// If user clicked the button, we want them to know we received their click
	if(manual) {
		busy = true;
		setTimeout(() => busy = false, 300);
	}

	// Validate criteria
	for(let i in criteria)
	{
		const criterion = criteria[i];
		try {
			for(let check of (criterion.checks || []))
				if(check.type == "file")
				{
					// Does file exist?
					if(check.action == "exists") {
						let stderr = null;
						await $CLI.exec(`ls ${check.path}`, d => stderr = d);
						if(stderr)
							throw "File not found";
						statuses[i] = true;
					}

					// Does file content match expectation? Define the right answer using a CLI invocation
					else if(check.action == "contents") {
						let stderr = null;
						await $CLI.exec(`ls ${check.path}`, d => stderr = d);
						if(stderr)
							throw "File not found";

						// Parse settings
						const commandExpected = check.commandExpected;
						const commandObserved = check.commandObserved || `cat ${check.path}`;

						// Calculate and compare expected vs. observed
						const expected = await $CLI.exec(commandExpected);
						const observed = await $CLI.exec(commandObserved);
						statuses[i] = observed == expected;
					}
				}
		} catch (error) {
			statuses[i] = false;
		}
	}
}
setTimeout(check, 500);
</script>
<div class="d-flex justify-content-between mt-4 mb-2">
	<div>
		<strong>Exercise Criteria:</strong>
	</div>
	<div>
		{statuses.filter(d => d).length} / {statuses.length}
	</div>
</div>

<ul class="list-group">
	{#each criteria as criterion, i}
		<li class="list-group-item list-group-item-action" class:list-group-item-success={statuses[i] === true}>
			<Icon name={statuses[i] ? "check-circle-fill" : "circle"} /> {@html criterion.name}
		</li>
	{/each}
</ul>

{#if hints.length > 0}
	<FormGroup class="mt-3" style="margin-bottom:0 !important"> <!-- FormGroup seems to add mb-3 automatically? -->
		<Label for="exampleRange">Showing {nbHints} out of {hints.length} hints</Label>
		<Input type="range" min={0} max={hints.length} step={1} bind:value={nbHints} />
	</FormGroup>
	<ListGroup class="mt-0">
		{#each hints.slice(0, nbHints) as hint}
			<ListGroupItem class="small">{@html hint}</ListGroupItem>
		{/each}
	</ListGroup>
{/if}

<button class="btn btn-sm btn-primary mt-3" class:btn-outline-primary={isDone} on:click={() => check(true)} disabled={isDone}>
	Check my work
	{#if busy}
		<Spinner size="sm" color="light" class="ms-2" />
	{/if}
</button>
