<script>
import { onMount } from "svelte";
import localforage from "localforage";
import { Button, Input } from "sveltestrap";
import Alert from "components/Alert.svelte";
import { tutorial } from "stores/tutorial";
import { getLocalForageKey } from "stores/config";

export let id;
export let choices = [];

// State
let checked = {};
let radio;
let error;
let success;

// State updates
$: stateId = getLocalForageKey("quiz") + `${$tutorial.id}-${$tutorial.step}-${id}`;
$: multiple = choices.filter(c => c.valid).length > 1;
$: if(id && (radio || checked)) updateState();

// Validate user response
function validate() {
	for(let i = 0; i < choices.length; i++) {
		const choice = choices[i].valid;
		const user = multiple ? checked[i] : radio === choices[i].value;
		if(!choice && user || choice && !user) {
			error = "That doesn't look right";
			return;
		}
	}

	success = true;
}

// Update state
async function updateState() {
	await localforage.setItem(stateId, multiple ? checked : radio);
}

// Set state of quiz
onMount(async () => {
	const state = await localforage.getItem(stateId);
	if(!state) return;
	if(multiple)
		checked = state;
	else
		radio = state;
	validate();
})
</script>

<Alert color={success ? "success" : "primary"}>
	<h6>
		<slot name="prompt">
			Prompt
		</slot>
	</h6>

	{#each choices as choice, i}
		{#if multiple}
			<Input type="checkbox" bind:checked={checked[i]} label={choice.value} disabled={success} />
		{:else}
			<Input type="radio" bind:group={radio} value={choice.value} label={choice.value} disabled={success} />
		{/if}
	{/each}

	<Button size="sm" color="info" class="mt-3" on:click={validate} disabled={success}>
		Submit
	</Button>

	{#if success}
		<p class="mt-2 text-success small">
			That is correct!
		</p>	
	{:else if error}
		<p class="mt-2 text-danger small">
			{error}
		</p>
	{/if}
</Alert>
