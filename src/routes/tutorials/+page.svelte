<script>
import { tutorials } from "$stores/tutorials";
import Listings from "$components/Listings.svelte";
</script>

<svelte:head>
	<title>Tutorials - sandbox.bio</title>
</svelte:head>

<Listings items={$tutorials} />
