<script>
import Execute from "$components/Execute.svelte";
</script>

Indexing a genome sorted BAM file allows one to quickly extract alignments
overlapping particular genomic regions. Moreover, indexing is required by
genome viewers such as IGV so that the viewers can quickly display alignments in each genomic region to which you navigate.

<Execute command="samtools index sample.sorted.bam" />

This will create an additional "index" file. List (`ls`) the contents of the current directory and look for the new index file: <Execute command="ls" inline />

Now, let's exploit the index to extract alignments from chromosome 20. To do
this, we use the samtools `view` command, which we will give proper treatment
in the next section. For now, just do it without understanding. No really. Do it.

<Execute command="samtools view sample.sorted.bam 20:1.4e6-1.5e6 | head -n 5" />

How many alignments are there in this region?

<Execute command="samtools view -c sample.sorted.bam 20:1.4e6-1.5e6" />
