<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

We've generated some paired random nucleotide data as gzip FASTQ using <Link href="https://github.com/noamteyssier/nucgen">`nucgen`</Link>.
These will be used to demonstrate the functionality of `bqtools`.

Let's take a look at what files we now have:

Type <Execute command="ls" inline /> in the command line.

You'll see a directory `./fastq` - let's take a look at the contents of this directory:

<Execute command="ls -lah fastq/" inline />

This directory contains 8 files, representing 4 different samples each with a pair of R1 and R2 files.

To get a sense of what these look like, let's explore the first sample:

<Execute command="zcat ./fastq/sample1_R1.fastq.gz | head -n 8" inline />

You'll notice that each record is a 4-line block of text:
1. sequence header
2. sequence
3. quality header (usually just a `+`)
4. quality scores
