<script>
import { Alert, Button, FormGroup, Input } from "sveltestrap";
import { supabaseAnon } from "$src/utils";

// Form state
let email = "";
let password = "";
let error;
let busy;

// Log in or sign up
async function login() {
	busy = true;
	const { error: err } = await supabaseAnon.auth.signInWithPassword({ email, password });
	error = err?.message;

	if (err) {
		busy = false;
	} else {
		error = false;
	}
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
			if (e.key === "Enter") login();
		}}
		disabled={busy}
	/>
</FormGroup>

{#if error}
	<Alert color="danger">
		{error}
	</Alert>
{/if}

<Button color="primary" class="w-100" on:click={login} disabled={busy}>Log in</Button>
