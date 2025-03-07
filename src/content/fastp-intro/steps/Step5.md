<script>
import Execute from "$components/Execute.svelte";
</script>

So far, we focused on using `fastp` to evaluate data quality. But if we want to use the filtered data in downstream analyses, we need to save the filtered FASTQ files for later processing. To do so, we can use the `--out` parameters:

<Execute command="fastp \ --in1 HG004_R1.fastq.gz \ --in2 HG004_R2.fastq.gz \ --out1 HG004_R1_filtered.fastq.gz \ --out2 HG004_R2_filtered.fastq.gz" />

This creates the original HTML/JSON reports, along with a filtered FASTQ file for each input: <Execute inline command="ls -lh" />. Notice that FASTQ files with `filtered` in their filenames are smaller, as expected.
