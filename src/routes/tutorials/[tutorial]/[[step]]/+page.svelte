<script>
import { page } from "$app/stores";
import { tutorial } from "$stores/tutorial";
import Tutorial from "$content/Tutorial.svelte";

$: id = $page.params.tutorial;
$: step = +$page.params.step || 0;
</script>

<svelte:head>
	{#if $tutorial.name}
		<title>{$tutorial.name} - sandbox.bio</title>
	{:else}
		<title>sandbox.bio</title>
	{/if}
</svelte:head>

{#if id}
	<Tutorial {id} {step} />
{/if}
