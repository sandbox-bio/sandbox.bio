<script>
import { CoreUtils } from "../terminal/coreutils";
import { Icon, Spinner } from "sveltestrap";

export let criteria = [];  // List of criteria that must be true for the exercise to be complete
let status = [];           // Status of each criteria (true/false)
let busy = false;       // Whether the user manually asked to check their work

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
		for(let check of (criterion.checks || []))
		{
			if(check.type == "file")
			{
				// Does file exist?
				if(check.action == "exists")
					status[i] = await CoreUtils.CLI.ls(check.path) !== false;
				// Does file content match expectation?
				else if(check.action == "contents")
				{
					const observed = await CoreUtils.CLI.cat(check.path);
					let expected;

					// If we define the right answer with a function we call
					if(check.fn)
						expected = await check.fn();
					// If we define the right answer using a CLI invocation
					else if(check.command)
						expected = await CoreUtils.CLI.exec(check.command);

					// Is it correct?
					if(check.output)
						await CoreUtils.FS.writeFile(check.output, expected);
					status[i] = observed == expected;
				}
			}
		}
	}

	// Check again
	setTimer();
}

// Check answers regularly
function setTimer() {
	setTimeout(check, 1000);
}

setTimer();
</script>
<p class="mt-4">
	<strong>Solution Criteria:</strong>
</p>

<ul class="list-group">
	{#each criteria as criterion, i}
		<li class="list-group-item list-group-item-action" class:list-group-item-success={status[i] === true}>
			<Icon name={status[i] ? "check-circle-fill" : "circle"} /> {@html criterion.name}
		</li>
	{/each}
</ul>

<button class="btn btn-sm btn-primary mt-3" on:click={() => check(true)} disabled={status.filter(d => d).length == status.length}>
	Check my work
	{#if busy}
		<Spinner size="sm" color="light" class="ms-2" />
	{/if}
</button>
