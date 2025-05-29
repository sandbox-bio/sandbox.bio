<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

Next, run:

<Execute command={"bowtie2 \ -x $REF \ -U reads_1.fq \ -S eg1.sam"} />

This runs the Bowtie 2 aligner, which aligns a set of unpaired reads to the <Link href="http://en.wikipedia.org/wiki/Lambda_phage">Lambda phage</Link> reference genome using the index generated in the previous step. The alignment results in SAM format are written to the file `eg1.sam`, and a short alignment summary is written to the console. (Actually, the summary is written to the "standard error" or "stderr" filehandle, which is typically printed to the console.)

To see the first few lines of the SAM output, run:

<Execute command="head -n 5 eg1.sam" />

The first few lines (beginning with `@`) are SAM header lines, and the rest of the lines are SAM alignments, one line per read or mate. See the <Link href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output">Bowtie 2 manual section on SAM output</Link> and the <Link href="http://samtools.sourceforge.net/SAM1.pdf">SAM specification</Link> for details about how to interpret the SAM file format.
