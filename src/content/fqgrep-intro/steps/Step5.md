<script>
import Execute from "$components/Execute.svelte";
</script>

Sometimes you want to find reads that **don't** contain a pattern — for example, filtering out adapter-contaminated reads or excluding certain barcodes.

## The -v flag

Use `-v` to invert the match — fqgrep will output reads that do **NOT** match the pattern:

<Execute command={`fqgrep -c 'AGATCGGAAGAGC' reads.fastq`} />

<Execute command={`fqgrep -v -c 'AGATCGGAAGAGC' reads.fastq`} />

Notice that the counts add up to the total number of reads in the file!

## Filtering out contaminants

Let's remove all reads containing contaminant sequences:

<Execute command={`fqgrep -v -f contaminants.txt reads.fastq | head -20`} />

You can save the clean reads to a new file (conceptually):

```bash
fqgrep -v -f contaminants.txt reads.fastq > clean_reads.fastq
```

## Excluding specific barcodes

If you want reads that don't have any of the known barcodes:

<Execute command={`fqgrep -v -c -f barcodes.txt reads.fastq`} />

## Multiple exclusions

You can combine `-v` with multiple patterns using `-e` or `-f`. This excludes reads matching **any** of the patterns:

<Execute command={`fqgrep -v -e 'AGATCGGAAGAGC' -e 'AGATGTGTATAAGAGACAG' -c reads.fastq`} />

This is equivalent to:

<Execute command={`fqgrep -v -f contaminants.txt -c reads.fastq`} />

> **Important**: With `-v`, a read is output only if it matches **none** of the specified patterns. This is perfect for cleaning data before analysis.
