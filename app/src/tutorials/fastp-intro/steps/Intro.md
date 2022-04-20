<script>
import Link from "components/Link.svelte";
</script>

Picture the scene: your sequencing run just finished, and you have FASTQ files in hand. Before diving into data analysis, let's start by inspecting the data quality. In this tutorial, we explore how to use the command-line tool <Link href="https://github.com/OpenGene/fastp">fastp</Link> to that end.

Let's start with the basics: `fastp` is a tool that is often used for two main purposes:

1. Generating a QC report to help you evaluate the quality of your sequencing data; and
2. Filtering out low quality data, whether it's throwing away entire reads or trimming those reads.
