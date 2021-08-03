<script>
import Execute from "components/Execute.svelte";
</script>

When you align FASTQ files with all current sequence aligners, the
alignments produced are in random order with respect to their position
in the reference genome. In other words, the BAM file is in the order
that the sequences occurred in the input FASTQ files.

<Execute command={"samtools view sample.bam | head -n 5"} />

Doing anything meaningful such as calling variants or visualizing
alignments in IGV) requires that the BAM is further manipulated. It must be sorted such that the alignments occur in "genome order". That is, ordered positionally based upon their alignment coordinates on each chromosome.

<Execute command={"samtools sort sample.bam -o sample.sorted.bam"} />

Now let's check the order:

<Execute command={"samtools view sample.sorted.bam | head -n 5"} />

Notice anything different about the coordinates of the alignments?
