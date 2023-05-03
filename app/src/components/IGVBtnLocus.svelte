<script>
import { Button } from "sveltestrap";
import { status } from "./stores/status";

export let locus = "all";
</script>

<Button color="outline-primary" on:click={() => $status.igv = { locus }}>
    Go to <code>{locus}</code>
</Button>
