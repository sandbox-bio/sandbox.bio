<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

<Link href="https://github.com/arcinstitute/binseq">`BINSEQ`</Link> is a file format for storing nucleotide sequences as binary data.
It is designed to be efficient and compact, making it suitable for large-scale genomic data analysis.

It has some benefits that make it an ergonomic choice for working with genomic data:
- Excellent compression of data
- Efficient parallel processing with random-access to records
- Paired records available in a single file (no separate R1/R2)

Our goal in this tutorial is to demonstrate how to use `bqtools`, a command-line tool that provides a set of utilities for working with BINSEQ files in the spirit of Unix utilities and other bioinformatics tools like `samtools` and `bedtools`.

We will cover the basics of importing/exporting sequencing data into BINSEQ files as well as some useful and common operations that can be performed with `bqtools` such as counting records, concatenation, and sequence grep.

For some more information about BINSEQ, its features, specifications, and how it compares to other formats, check out the [BINSEQ preprint](https://www.biorxiv.org/content/10.1101/2025.04.08.647863v2).

For more complete information about `bqtools`, check out its [documentation](https://github.com/arcinstitute/bqtools)

For programmatic access to BINSEQ, check out the [BINSEQ rust docs](https://docs.rs/binseq/latest/binseq/).
