<script>
import Execute from "components/Execute.svelte";
import IGV from "components/IGV.svelte";

let isOpen = false;
let igvOptions = {
	locus: "chr20:1,299,889-1,300,567",
	tracks: [
		{ url: "./data/samtools-intro/sample.bam", name: "Read alignment" },
	]
};
</script>

We preloaded sample data into your sandbox as `sample.sam`. Type <Execute command={"ls"} inline={true} /> in the command line to see it.

To get a rough sense of what the alignments in `sample.sam` look like, let's load them into IGV: <button class="btn btn-sm btn-primary" on:click={() => isOpen = !isOpen}>Launch IGV</button>

<IGV options={igvOptions} bind:isOpen={isOpen} />
