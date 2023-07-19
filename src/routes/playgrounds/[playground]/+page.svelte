<script>
import { page } from "$app/stores";
import Sandbox from "$components/playgrounds/Sandbox.svelte";
import { TOOLS } from "$stores/sandbox";

$: playground = $page.params.playground;
$: tool = TOOLS.find((d) => d.name === playground);
</script>

{#if tool}
	<Sandbox {tool} />
{:else}
	Tool not recognized.
{/if}
