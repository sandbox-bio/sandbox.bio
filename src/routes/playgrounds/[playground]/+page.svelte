<script>
import { onMount } from "svelte";
import Aioli from "@biowasm/aioli";
import { page } from "$app/stores";
import Sandbox from "$components/playgrounds/Sandbox.svelte";
import { sandbox, TOOLS } from "$stores/sandbox";

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
	CLI = await new Aioli(config, { printInterleaved: false });
	ready = true;
});
</script>

Playgrounds

{#if tool}
	<Sandbox {tool} {CLI} />
{:else}
	<p>Playground for tool <code>{playground}</code> does not exist.</p>

	<p>Please <a href="https://github.com/sandbox-bio/sandbox.bio/discussions">reach out</a> if you would like to propose new playgrounds!</p>
{/if}
