<script>
import { onMount } from "svelte";
import { Button } from "sveltestrap";
import Aioli from "@biowasm/aioli";
import { page } from "$app/stores";
import { sandbox, TOOLS } from "$stores/sandbox";
import Sandbox from "$components/playgrounds/Sandbox.svelte";

// State
let CLI;
let ready = false; // false if still initializing Aioli

$: playground = $page.params.playground;
$: tool = TOOLS.find((d) => d.name === playground);

// Initialize on page load
onMount(async () => {
	// Create store with data from localforage
	await sandbox.init();

	// Initialize all CLI tools
	const config = TOOLS.map((d) => d.aioli);
	console.log("config", config);
	CLI = await new Aioli(
		[
			{ tool: "base", version: "1.0.0" }, // need at least 1 tool with reinit=false,
			...config
		],
		{ printInterleaved: false }
	);
	ready = true;
});
</script>

<svelte:head>
	<title>{tool.name} playground - sandbox.bio</title>
</svelte:head>

{#if tool}
	<Sandbox {tool} {CLI} {ready}>
		<div slot="playgrounds">
			{#each TOOLS as toolLink}
				<Button
					href="/playgrounds/{toolLink.name}"
					size="sm"
					class="ms-2"
					color={toolLink.name === tool.name ? "outline-secondary" : "outline-primary"}
				>
					{toolLink.name}
				</Button>
			{/each}
		</div>
	</Sandbox>
{:else}
	<p>Playground for tool <code>{playground}</code> does not exist.</p>

	<p>Please <a href="https://github.com/sandbox-bio/sandbox.bio/discussions">reach out</a> if you would like to propose new playgrounds!</p>
{/if}
