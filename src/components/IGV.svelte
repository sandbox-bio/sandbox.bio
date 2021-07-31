<script>
import { Spinner, Modal } from "sveltestrap";

// State
export let options = {};    // IGV.js options
export let isOpen = false;  // Whether modal is showing or not
let igvPromise;             // Resolves when igv.js is done loading

function createIGV() {
	var igvDiv = document.getElementById("igv-div");
	igvPromise = igv.createBrowser(igvDiv, options);
}
</script>

<svelte:head>
	<script src="https://cdn.jsdelivr.net/npm/igv@2.9.3/dist/igv.min.js"></script>
</svelte:head>

<Modal body header="IGV.js" size="xl" on:open={createIGV} toggle={() => isOpen = !isOpen} {isOpen}>
	{#await igvPromise}
		<Spinner size="sm" color="primary" type="border" /> Loading...
	{/await}

	<slot name="before"></slot>

	<div id="igv-div"></div>
	<br />
	<slot name="after"></slot>
</Modal>
