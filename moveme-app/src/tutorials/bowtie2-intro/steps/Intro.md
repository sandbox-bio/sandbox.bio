<script>
import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";
</script>

<Alert>
	**Note**: This tutorial is an interactive version of the <Link href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example">bowtie2 tutorial</Link> developed by the <Link href="https://langmead-lab.org/">Langmead Lab</Link>. The contents are the same, but the data was subsampled so it can be analyzed in your browser.
</Alert>

In this tutorial, we'll explore how to use `bowtie2` to align DNA sequencing reads to the Lambda reference genome.

We'll see examples of aligning both single-end and paired-end reads stored in `FASTQ` files.
