<script>
import { onMount } from "svelte";
import { status } from "./stores/status";

export let options = {};
let loading = true;
let igvDiv;
let browser = {};

onMount(async() => {
	// Explicitly specify the reference URLs so that igv.js doesn't try downloading RefSeq genes
	if(!options.genome) {
		options.reference = {
			id: "hg19",
			name: "Human (hg19)",
			fastaURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta",
			indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta.fai",
			cytobandURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/cytoBand.txt",
			tracks: []
		};
	}

	// Create IGV browser
	browser = await igv.createBrowser(igvDiv, options);
	loading = false;
});

// Allow tutorials to interactively change IGV status
$: if($status.igv) {
	// Note that `browser.currentLoci` can give fractional coordinates
	const locusCurrent = browser.referenceFrameList.map((locus) => locus.getLocusString()).join(" ");

	// Locus change
	if($status.igv.locus && locusCurrent !== $status.igv.locus) {
		browser.search($status.igv.locus);
	}

	// Add a new track
	if($status.igv.loadTrack) {
		browser.loadTrack($status.igv.loadTrack);
	}
}
</script>

<div style="opacity: { loading ? 0.6 : 1 }">
	<div bind:this={igvDiv} />
</div>
