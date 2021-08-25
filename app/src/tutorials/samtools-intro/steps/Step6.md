<script>
import Execute from "components/Execute.svelte";
import IGV from "components/IGV.svelte";

let isOpen = false;
let igvOptions = {
	locus: "chr20:1,299,889-1,300,567",
	tracks: [
		{ url: "https://storage.googleapis.com/sandbox.bio/data/sample.bam", name: "Read alignment" },
	]
};
</script>

The samtools `view` command is the most versatile tool in the samtools package.
It's main function, not surprisingly, is to allow you to convert the binary
(i.e., easy for the computer to read and process) alignments in the BAM file
view to text-based SAM alignments that are easy for *humans* to read and process.

Let us start by inspecting the first five alignments in our BAM in detail:

<Execute command={"samtools view sample.sorted.bam | head -n 5"} />

For each read, can you identify where in the genome the read landed? With what mapping quality? Can you parse what the SAM flags mean? (use <Execute command={"samtools flags"} inline={true} />)

Finally, let's visualize the alignments: <button class="btn btn-sm btn-primary" on:click={() => isOpen = !isOpen}>Launch IGV</button>

<IGV options={igvOptions} bind:isOpen={isOpen} />

Using IGV, can you visually identify reads where the DNA sequence differs from that of the reference?
