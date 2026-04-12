<script>
import Execute from "$components/Execute.svelte";
</script>

One of fqgrep's most powerful features is proper handling of paired-end reads. When one mate matches a pattern, you usually want **both** mates in the output to keep pairs synchronized.

## The problem with grep

If you use regular `grep` on paired-end files separately, you'll get different reads from R1 and R2 — breaking the pairing that downstream tools require.

## The --paired flag

fqgrep's `--paired` flag treats input files as paired. When a pattern matches in **either** R1 or R2, **both** reads are output (interleaved):

<Execute command={`fqgrep --paired 'AGATCGGAAGAGC' reads_R1.fastq reads_R2.fastq | head -24`} />

Notice the output alternates between R1 and R2 reads (interleaved), keeping pairs together.

## Counting paired matches

<Execute command={`fqgrep --paired -c 'AGATCGGAAGAGC' reads_R1.fastq reads_R2.fastq`} />

This counts **read pairs** where at least one mate contains the adapter.

## Why this matters

Our sample data has adapter contamination in different mates:

- Some pairs have adapter in R1 only
- Some pairs have adapter in R2 only

With `--paired`, fqgrep finds all contaminated pairs regardless of which mate has the adapter:

<Execute command={`fqgrep --paired --color always 'AGATCGGAAGAGC' reads_R1.fastq reads_R2.fastq | head -16`} />

## Searching for Tn5 in paired data

Let's find pairs where either mate contains the Tn5 mosaic end:

<Execute command={`fqgrep --paired -c 'AGATGTGTATAAGAGACAG' reads_R1.fastq reads_R2.fastq`} />

## A note on -v with --paired

Be careful combining `-v` with `--paired`. Remember that `--paired` outputs a pair if **either** mate is selected. With `-v`, a mate is selected if it does **not** match. So `--paired -v` outputs pairs where **at least one** mate doesn't match — which is almost all pairs, not the "clean pairs" you might expect.

To truly filter out contaminated pairs, you'd need to run fqgrep without `--paired` on each file separately and intersect the results.

> **Key insight**: The `--paired` flag is essential for maintaining read pair integrity. Without it, you'd need complex scripts to re-synchronize filtered R1 and R2 files.
