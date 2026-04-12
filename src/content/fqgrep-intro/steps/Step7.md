<script>
import Execute from "$components/Execute.svelte";
</script>

DNA is double-stranded, and sequencing can capture either strand. A sequence you're looking for might appear in its reverse complement form depending on the read orientation.

## The --reverse-complement flag

fqgrep can automatically search both the original pattern AND its reverse complement using `--reverse-complement`.

Let's see the difference with the Illumina adapter. First, a forward-only search:

<Execute command={`fqgrep -c 'AGATCGGAAGAGC' reads.fastq`} />

Now searching both orientations:

<Execute command={`fqgrep --reverse-complement -c 'AGATCGGAAGAGC' reads.fastq`} />

The reverse complement search finds more reads because some contain the adapter in the opposite orientation (`GCTCTTCCGATCT`). Without `--reverse-complement`, those reads would be missed.

## Why this matters

Consider an adapter or primer sequence. Depending on:

- Which strand was sequenced
- Read 1 vs Read 2 orientation
- Library prep protocol

...the sequence might appear as its reverse complement. Searching both orientations ensures you find all relevant reads.

## Visual inspection

Let's see matches with color highlighting — notice some reads match the forward adapter and others match the reverse complement:

<Execute command={`fqgrep --reverse-complement --color always 'AGATCGGAAGAGC' reads.fastq | head -16`} />

## A special case: palindromic sequences

Some sequences are their own reverse complement. For example, `ACGTACGT` is a palindrome:

<Execute command={`fqgrep -c 'ACGTACGT' reads.fastq`} />

<Execute command={`fqgrep --reverse-complement -c 'ACGTACGT' reads.fastq`} />

The counts are the same because the forward and reverse complement sequences are identical — `--reverse-complement` has no effect for palindromes.

## Combining with paired-end

You can use `--reverse-complement` with `--paired` for comprehensive paired-end searches:

<Execute command={`fqgrep --paired --reverse-complement -c 'AGATCGGAAGAGC' reads_R1.fastq reads_R2.fastq`} />

This finds pairs where either mate contains the adapter in either orientation.

> **Pro tip**: When searching for adapter or primer sequences, always consider whether you need `--reverse-complement`. For adapters that ligate to both ends of inserts, searching both orientations is usually appropriate.
