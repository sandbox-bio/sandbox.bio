<script>
import { createEventDispatcher } from "svelte";
import { Alert } from "sveltestrap";

export let type = "login";  // Supported: login, signup
export let error = false;
export let success = false;

// Form state
let email = "";
let password = "";
const dispatch = createEventDispatcher();

// Log in or sign up
function run() {
	dispatch(type, {
		email, password
	});
}
</script>

<div class="form-floating mb-1">
	<input type="email" class="form-control rounded-4" id="email-{type}" placeholder="name@example.com" bind:value={email}>
	<label for="email-{type}">Email address</label>
</div>

<div class="form-floating mb-3">
	<input type="password" class="form-control rounded-4" id="password-{type}" placeholder="Password" bind:value={password} on:keypress={e => { if (e.key === "Enter")run(); }}>
	<label for="password-{type}">Password</label>
</div>

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

<button class="w-100 mb-2 btn btn-lg rounded-4 btn-primary" on:click={run}>
	{#if type == "login"}
		Log in
	{:else}
		Create account
	{/if}
</button>
