<script>
import { tutorials } from "$stores/tutorials";
import Listings from "$components/Listings.svelte";
</script>

<Listings items={$tutorials} />
