<script>
import Execute from "$components/Execute.svelte";
</script>

Often you need to search for several patterns at once — multiple barcodes, different adapter variants, or a set of sequences of interest.

## Using -e for multiple patterns

The `-e` (or `--regexp`) option lets you specify multiple patterns. A read matches if it contains **any** of the patterns:

<Execute command={`fqgrep -e 'ACGTACGT' -e 'TGCATGCA' -e 'GCTAGCTA' reads.fastq | head -20`} />

Let's count how many reads match any of our three barcodes:

<Execute command={`fqgrep -c -e 'ACGTACGT' -e 'TGCATGCA' -e 'GCTAGCTA' reads.fastq`} />

## Using a pattern file with -f

When you have many patterns, it's easier to put them in a file. Check out our barcodes file:

<Execute command={`cat barcodes.txt`} />

Now search using `-f` (or `--file`):

<Execute command={`fqgrep -f barcodes.txt reads.fastq | head -20`} />

<Execute command={`fqgrep -c -f barcodes.txt reads.fastq`} />

## Contaminant screening

A common use case is screening for multiple contaminant sequences. Our `contaminants.txt` file contains both the Illumina adapter and Tn5 mosaic end:

<Execute command={`cat contaminants.txt`} />

<Execute command={`fqgrep -c -f contaminants.txt reads.fastq`} />

This gives you a quick count of reads with any contaminant sequence.

## Combining -e and -f

You can even combine command-line patterns with a pattern file:

<Execute command={`fqgrep -f barcodes.txt -e 'AGATCGGAAGAGC' -c reads.fastq`} />

This searches for all barcodes AND the adapter pattern.

> **Pro tip**: For quality control, create a standard contaminants file for your lab/pipeline and use it to quickly screen new sequencing runs.
