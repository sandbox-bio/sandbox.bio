<script>
import { createEventDispatcher } from "svelte";
import { Alert, Button, FormGroup, Input } from "sveltestrap";

export let type = "login"; // Supported: login, signup
export let error = false;
export let success = false;
export let busy = false;

// Form state
let email = "";
let password = "";
const dispatch = createEventDispatcher();

// Log in or sign up
function run() {
	dispatch(type, {
		email,
		password
	});
}
</script>

<FormGroup floating label="Email address">
	<Input type="email" bind:value={email} disabled={busy} />
</FormGroup>

<FormGroup floating label="Password">
	<Input
		type="password"
		bind:value={password}
		on:keypress={(e) => {
			if (e.key === "Enter") run();
		}}
		disabled={busy}
	/>
</FormGroup>

{#if error}
	<Alert color="danger">
		{error}
	</Alert>
{/if}

{#if success}
	<Alert color="success">
		{success}
	</Alert>
{/if}

<Button color="primary" class="w-100" on:click={run} disabled={busy}>
	{#if type == "login"}
		Log in
	{:else}
		Create account
	{/if}
</Button>
