<script>
import { status } from "status";
import { CLI } from "terminal/cli";
import { Icon, Spinner } from "sveltestrap";

export let criteria = [];        // List of criteria that must be true for the exercise to be complete
let statuses = [];                 // Status of each criteria (true/false)
let busy = false;                // Whether the user manually asked to check their work

// Check status every time a command finishes running
$: if($status.terminal == "execDone"){
	console.log("Checking solutions...")
	check();
};

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
						await $CLI.exec(`ls ${check.path}`);
						statuses[i] = true;
					}

					// Does file content match expectation? Define the right answer using a CLI invocation
					else if(check.action == "contents") {
						const observed = await $CLI.exec(`cat ${check.path}`);
						const expected = await $CLI.exec(check.command);
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
<p class="mt-4">
	<strong>Solution Criteria:</strong>
</p>

<ul class="list-group">
	{#each criteria as criterion, i}
		<li class="list-group-item list-group-item-action" class:list-group-item-success={statuses[i] === true}>
			<Icon name={statuses[i] ? "check-circle-fill" : "circle"} /> {@html criterion.name}
		</li>
	{/each}
</ul>

<button class="btn btn-sm btn-primary mt-3" on:click={() => check(true)} disabled={statuses.filter(d => d).length == statuses.length && statuses.length != 0}>
	Check my work
	{#if busy}
		<Spinner size="sm" color="light" class="ms-2" />
	{/if}
</button>
