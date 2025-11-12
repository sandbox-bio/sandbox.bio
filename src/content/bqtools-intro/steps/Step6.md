<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

Now we have a few BINSEQ files - what can we do with them?

If we want to explore the underlying sequence data we can make use of the `binseq decode` command.
This will let us output the sequence data in a human-readable format.

Let's output a few sequences as a tab-delimited file (TSV).

<Execute command="bqtools decode ./fastq/sample1.vbq | head -n 20" />

Notice that the first column is the record ID and the second column is the sequence itself.
You'll notice that the record ID is repeated every two lines for this file - this is because the R1 and R2 sequence data currently have the same record ID (its position in the file).

> Note: You may notice that if you run this command multiple times, the record IDs positions may change. BINSEQ files are meant to be processed in parallel, and as such they are not guaranteed to be output in the same order every time. You may process them sequentially by setting the `-T/--threads` flag to 1.

We can also decode these files to FASTA, FASTQ, and even BAM by making use of the `-f/--format` flag.

Let's take a look at a FASTA representation of the sequences.

<Execute command="bqtools decode ./fastq/sample1.vbq -fa | head -n 20" />

And again as FASTQ:

<Execute command="bqtools decode ./fastq/sample1.vbq -fq | head -n 20" />
