<script>
import { page } from "$app/stores";
import Tutorial from "$content/Tutorial.svelte";

$: tutorial = $page.params.tutorial;
$: step = $page.params.step || 0;
</script>

{#if tutorial}
	<Tutorial id={tutorial} {step} />
{/if}
