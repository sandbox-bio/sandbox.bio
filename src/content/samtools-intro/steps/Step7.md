<script>
import Execute from "$components/Execute.svelte";
</script>

The "header" in a BAM file records important information regarding the 
reference genome to which the reads were aligned, as well as other information
about how the BAM has been processed. One can ask the `view`
command to report solely the header by using the `-H` option.

<Execute command={"samtools view -H sample.sorted.bam"} />
