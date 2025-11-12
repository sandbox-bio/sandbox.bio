<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

<Link href="https://github.com/arcinstitute/binseq">`BINSEQ`</Link> is a file format for storing nucleotide sequences as binary data.
It is designed to be efficient and compact, making it suitable for large-scale genomic data analysis.

To convert a FASTQ (or FASTA) file to BINSEQ we can make use of the `bqtools` command-line tool.
Specifically we can use `bqtools encode` to convert the file.

Lets encode a single FASTQ file:

<Execute command="bqtools encode fastq/sample1_R1.fastq.gz -o sample1_R1.vbq"/>

This will run a BINSEQ conversion using default parameters and create a new `.vbq` file in the current directory.

Let's compare the file sizes of the original FASTQ file and the resulting BINSEQ file:

<Execute command="ls -lh fastq/sample1_R1.fastq.gz sample1_R1.vbq"/>

You'll see that that the BINSEQ file is significantly smaller than the original FASTQ file!

By default, `bqtools encode` will write a VBQ file.
This is one of the two variants of BINSEQ and is the more flexible variants that supports variable-length sequences, quality scores, and sequence headers - making it fully lossless to FASTQ.

We can also write a BQ file - the simpler BINSEQ variant that only supports fixed-length sequences without quality scores or headers.
This format is more limited than the VBQ format, but can lead to massive throughput improvements for applications that do not require quality scores or headers.

`bqtools` can auto-determine the mode written using the provided extension (`.vbq` or `.bq`)

<Execute command="bqtools encode fastq/sample1_R1.fastq.gz -o sample1_R1.bq"/>

> Note: you can also set the BINSEQ mode (bq/vbq) using the `-m/--mode` option:


Let's compare the file sizes of all the files:

<Execute command="ls -lh fastq/sample1_R1.fastq.gz sample1_R1.vbq sample1_R1.bq"/>

> Note: In this case, on very small files, the BQ file is slightly smaller than the VBQ file.

Because BQ does not use any compression (just two-bit encoding of nucleotides), it has a deterministic size from the number of nucleotides in the file.

VBQ uses <Link href="https://en.wikipedia.org/wiki/Zstd">ZSTD</Link> encoding internally and depending on the nucleotide characteristics, the file size is non-deterministic but tends to be smaller than BQ in most cases and roughly the same size as <Link href="https://en.wikipedia.org/wiki/CRAM_(file_format)">CRAM</Link> files.
