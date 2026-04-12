<script>
import Execute from "$components/Execute.svelte";
</script>

**fqgrep** is a fast, specialized tool for searching patterns in FASTQ files. Developed by [Fulcrum Genomics](https://www.fulcrumgenomics.com/), it's like `grep` but designed specifically for sequencing data.

## Why not just use grep?

You might think: "Why can't I just use regular `grep` to search my FASTQ files?" Here's the problem:

FASTQ files have a specific 4-line structure for each read:

1. Header line (starts with `@`)
2. Sequence line
3. Plus line (`+`)
4. Quality line

When you use `grep` to find a sequence pattern, you only get the matching line — not the complete FASTQ record. Worse, using `grep -B1 -A2` to grab surrounding lines is extremely slow and error-prone.

## fqgrep advantages

**fqgrep** solves these problems:

- **FASTQ-aware**: Always outputs complete, valid FASTQ records
- **Blazing fast**: Up to ~1000x faster than `grep` ([benchmarks](https://github.com/Rbfinch/grepq))
- **Paired-end support**: Keeps read pairs together when one mate matches
- **grep-compatible**: Familiar options like `-e`, `-f`, `-v`, `-c`, `-F`

## Getting started

Let's check that fqgrep is available:

<Execute command={`fqgrep --version`} />

Take a look at the available options:

<Execute command={`fqgrep --help`} />

## Sample data

This tutorial includes synthetic FASTQ files with various patterns planted in them:

- **reads.fastq** — Single-end reads with adapters, barcodes, and other sequences
- **reads_R1.fastq** and **reads_R2.fastq** — Paired-end reads

Let's see what we're working with:

<Execute command={`head -12 reads.fastq`} />

This shows the first 3 reads (4 lines each). Notice the structure: header, sequence, `+`, quality scores.

> fqgrep understands FASTQ structure, so it always returns complete records — never partial reads that would break downstream tools.
