<script>
import { Modal } from "sveltestrap";
import IGV from "./IGV.svelte";

// State
export let options = {}; // IGV.js options
export let isOpen = false; // Whether modal is showing or not
</script>

<Modal body header="IGV.js" size="xl" toggle={() => (isOpen = !isOpen)} {isOpen}>
	<IGV {options} />

	<slot name="after" />
</Modal>
