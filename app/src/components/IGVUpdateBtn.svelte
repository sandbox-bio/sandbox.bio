<script>
import { Button } from "sveltestrap";
import { status } from "./stores/status";

export let locus = null;
export let loadTrack = null;
</script>

<Button class="mb-3" color="primary" on:click={() => $status.igv = { locus, loadTrack }}>
	<slot />
</Button>
