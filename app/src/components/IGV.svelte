<script>
import { onMount } from "svelte";
import { Spinner } from "sveltestrap";

export let options = {};
let loading = true;
let igvDiv;

onMount(async() => {
	await igv.createBrowser(igvDiv, options);
	loading = false;
});
</script>

{#if loading}
	<Spinner size="sm" color="primary" type="border" /> Loading...
{/if}

<div bind:this={igvDiv} />
