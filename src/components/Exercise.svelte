<script>
import { onDestroy, onMount } from "svelte";
import { Button, FormGroup, Icon, ListGroup, ListGroupItem, Spinner } from "sveltestrap";
import { EXEC_MODE_BUS, cli } from "$stores/cli";
import { DIR_TUTORIAL } from "$src/config";

export let criteria = []; // List of criteria that must be true for the exercise to be complete
export let hints = []; // Hints to show user

let statuses = []; // Status of each criteria (true/false)
let busy = false; // Whether the user manually asked to check their work
let working = false; // Whether we're currently checking exercise status (don't do it more than once at a time)
let nbHints = 0; // How many hints we're displaying
let timers = []; // Track setTimeout timers so we can stop them when done

$: isDone = statuses.filter((d) => d).length == statuses.length && statuses.length != 0;

onMount(() => {
	// Need to initialize `statuses` because `isDone` will be checked before all statuses are done
	// because we can't make `$cli.exec` a synchronous call.
	statuses = new Array(criteria.length).fill(false);
	timers.push(setTimeout(check, 500));
});

// When move away, stop all exercise timers
onDestroy(() => {
	timers.forEach((timer) => clearTimeout(timer));
});

// Validate user's input
async function check(manual = false) {
	// If user clicked the button, we want them to know we received their click
	if (manual) {
		busy = true;
		timers.push(setTimeout(() => (busy = false), 100));
	}

	// Validate criteria
	if (!working) {
		working = true;
		try {
			for (let i in criteria) {
				const criterion = criteria[i];
				try {
					// Does file exist?
					for (let check of criterion.checks || []) {
						if (check.action == "exists") {
							const ls = $cli.ls(`${DIR_TUTORIAL}/${check.path}`);
							if (!ls) throw "File not found";
							statuses[i] = true;
						}

						// Does file content match expectation? Define the right answer using a CLI invocation
						else if (check.action == "contents") {
							const ls = $cli.ls(`${DIR_TUTORIAL}/${check.path}`);
							if (!ls) throw "File not found";

							// Parse settings
							const commandExpected = check.commandExpected;
							const commandObserved = check.commandObserved || `cat ${check.path}`;

							// Diff expected vs. observed (diff outputs nothing if equal)
							await $cli.clearCache();
							$cli.exec(`cd ${DIR_TUTORIAL} && diff -q <(${commandObserved}) <(${commandExpected}) | wc -l`, {
								mode: EXEC_MODE_BUS,
								callbackExercise: (s) => {
									console.log("[Exercise check] Valid =", s == 0, s);
									statuses[i] = s == 0;
								}
							});
						}
					}
				} catch (error) {
					statuses[i] = false;
				}
			}
		} catch (error) {
			console.error(error);
		}

		// Check exercise status regularly
		if (!isDone) {
			// Manual check is one-time only
			if (!manual) timers.push(setTimeout(check, 1000));
		} else {
			console.warn(`Stopped checking exercises because exercise is complete.`);
		}

		working = false;
	}
}
</script>

<div class="d-flex justify-content-between mt-4 mb-2">
	<div>
		<strong>Exercise Criteria:</strong>
	</div>
	<div>
		{statuses.filter((d) => d).length} / {statuses.length}
	</div>
</div>

<ul class="list-group">
	{#each criteria as criterion, i}
		<li class="list-group-item list-group-item-action" class:list-group-item-success={statuses[i] === true}>
			<Icon name={statuses[i] ? "check-circle-fill" : "circle"} />
			{@html criterion.name}
		</li>
	{/each}
</ul>

<button class="mt-2 mb-4 btn btn-sm btn-primary" on:click={() => check(true)} disabled={isDone}>
	Check my work
	{#if busy}
		<Spinner size="sm" color="light" class="ms-2" />
	{/if}
</button>

{#if hints.length > 0}
	<FormGroup class="mt-3" style="margin-bottom:10px !important">
		<!-- FormGroup seems to add mb-3 automatically? -->
		<strong>Hints:</strong><br />
		<span class="small text-muted">Showing {nbHints}/{hints.length} hints</span>
		<Button size="sm" color="outline-primary" on:click={() => nbHints++} disabled={nbHints === hints.length}>Show more</Button>
		<Button size="sm" color="outline-primary" on:click={() => nbHints--} disabled={nbHints === 0}>Show fewer</Button>
	</FormGroup>
	<ListGroup class="mt-0">
		{#each hints as hint, i}
			<ListGroupItem class="small">
				{#if i < nbHints}
					{@html hint}
				{:else}
					<span class="text-muted"> [Hint hidden] </span>
				{/if}
			</ListGroupItem>
		{/each}
	</ListGroup>
{/if}
