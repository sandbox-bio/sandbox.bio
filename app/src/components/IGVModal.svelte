<script>
import { Modal } from "sveltestrap";
import IGV from "./IGV.svelte";

// State
export let options = {};    // IGV.js options
export let isOpen = false;  // Whether modal is showing or not

// Explicitly specify the reference URLs so that igv.js doesn't try downloading RefSeq genes
options.reference = {
	id: "hg19",
	name: "Human (hg19)",
	fastaURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta",
	indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta.fai",
	cytobandURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/cytoBand.txt",
	tracks: []
};
</script>

<Modal body header="IGV.js" size="xl" toggle={() => isOpen = !isOpen} {isOpen}>
	<IGV {options} />

	<slot name="after" />
</Modal>
