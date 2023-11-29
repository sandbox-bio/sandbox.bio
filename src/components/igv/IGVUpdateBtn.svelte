<script>
import { Button } from "sveltestrap";
import { status } from "$stores/status";

export let locus = null;
export let loadTrack = null;
</script>

<div class="text-center">
	<Button class="mb-3" color="primary" size="sm" on:click={() => ($status.igv = { locus, loadTrack })}>
		<slot>
			Go to {locus}
		</slot>
	</Button>
</div>
