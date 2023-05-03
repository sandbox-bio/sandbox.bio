<script>
import { onMount } from "svelte";
import { Spinner } from "sveltestrap";

export let options = {};
let loading = true;
let igvDiv;

onMount(async() => {
	// Explicitly specify the reference URLs so that igv.js doesn't try downloading RefSeq genes
	options.reference = {
		id: "hg19",
		name: "Human (hg19)",
		fastaURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta",
		indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta.fai",
		cytobandURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/cytoBand.txt",
		tracks: []
	};

	// Create IGV browser
	await igv.createBrowser(igvDiv, options);
	loading = false;
});
</script>

{#if loading}
	<Spinner size="sm" color="primary" type="border" /> Loading...
{/if}

<div bind:this={igvDiv} />
