<script>
import { Button, Icon, Input } from "sveltestrap";

let email;
let done = false;
let error;

async function subscribe() {
	error = null;
	const response = await fetch("/api/v1/mailinglist", {
		method: "POST",
		body: JSON.stringify({ email })
	}).then((d) => d.json());
	error = response.error;
	if (!error) {
		done = true;
	}
}
</script>

<!-- Using <Badge> looks weird on mobile -->
<div class="text-center d-flex">
	<div class="mx-auto bg-primary-subtle rounded-3 py-3 px-md-3 px-2 mt-4" style="font-size: 17px">
		<Icon name="bell-fill" />
		Get notified when new tutorials are available on sandbox.bio:
		<Input
			type="email"
			bind:value={email}
			disabled={done}
			class="mt-1"
			placeholder="Your email"
			on:keydown={(e) => {
				if (e.key === "Enter") {
					subscribe();
				}
			}}
		/>
		<Button color="primary" size="sm" class="mt-2" disabled={done} on:click={subscribe}>Subscribe</Button>
		{#if error}
			<p class="mt-2 mb-0 text-danger small">{error}</p>
		{/if}

		{#if done}
			<p class="mt-2 mb-0 text-success small">You've been subscribed!</p>
		{/if}
	</div>
</div>
