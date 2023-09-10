<script>
import Listings from "$components/Listings.svelte";
import { playgrounds } from "$stores/tutorials";
</script>

<Listings colLg={6} colXxl={3} items={$playgrounds} title="Playgrounds" />
