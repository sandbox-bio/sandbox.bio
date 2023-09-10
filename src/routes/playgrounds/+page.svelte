<script>
import Listings from "$components/Listings.svelte";
import { playgrounds } from "$stores/tutorials";
</script>

<Listings items={$playgrounds} title="Playgrounds" />
