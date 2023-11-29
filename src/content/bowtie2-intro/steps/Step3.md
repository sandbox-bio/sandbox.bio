<script>
import Execute from "$components/Execute.svelte";
</script>

To align paired-end reads included with Bowtie 2, stay in the same directory and run:

<Execute command="bowtie2 \ -x $REF \ -1 reads_1.fq \ -2 reads_2.fq \ -S eg2.sam" />

This aligns a set of paired-end reads to the reference genome, with results written to the file `eg2.sam`.
