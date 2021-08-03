<script>
import { Spinner, Modal } from "sveltestrap";

// State
export let options = {};    // IGV.js options
export let isOpen = false;  // Whether modal is showing or not
let igvPromise;             // Resolves when igv.js is done loading

function createIGV() {
	var igvDiv = document.getElementById("igv-div");

	// Explicitly specify the reference URLs so that igv.js doesn't try downloading RefSeq genes
	options.reference = {
		id: "hg19",
		name: "Human (hg19)",
		fastaURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta",
		indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta.fai",
		cytobandURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/cytoBand.txt",
		tracks: []
	};
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
