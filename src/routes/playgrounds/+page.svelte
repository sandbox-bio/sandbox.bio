<script>
import Listings from "$components/Listings.svelte";
import { playgrounds } from "$stores/tutorials";
</script>

<svelte:head>
	<title>Playgrounds - sandbox.bio</title>
</svelte:head>

<Listings items={$playgrounds} title="Playgrounds" />
