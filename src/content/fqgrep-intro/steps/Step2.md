<script>
import Execute from "$components/Execute.svelte";
</script>

The simplest use of fqgrep is searching for a specific sequence in your reads.

## Finding a pattern

Let's search for the Illumina TruSeq adapter sequence `AGATCGGAAGAGC`, which commonly appears at the end of reads when inserts are shorter than the read length:

<Execute command={`fqgrep 'AGATCGGAAGAGC' reads.fastq | head -20`} />

Notice that fqgrep returns **complete FASTQ records** — all 4 lines for each matching read.

## Counting matches

Often you just want to know _how many_ reads contain a pattern, not see all of them. Use `-c` (or `--count`):

<Execute command={`fqgrep -c 'AGATCGGAAGAGC' reads.fastq`} />

This is much faster than piping to `wc -l` because fqgrep doesn't need to format the output.

## Fixed strings vs. regex

By default, fqgrep treats patterns as **regular expressions**. This means characters like `.`, `*`, `+`, and `[` have special meaning. Use `-F` (or `--fixed-strings`) when you want to match the pattern **literally**.

Consider the pattern `GACG.GATTA`. In regex mode, the `.` matches **any** character, so this finds both `GACGAGATTA` and `GACGTGATTA`:

<Execute command={`fqgrep -c 'GACG.GATTA' reads.fastq`} />

With `-F`, fqgrep validates that fixed-string patterns contain only valid DNA bases — so non-DNA characters like `.` are rejected:

<Execute command={`fqgrep -F -c 'GACGAGATTA' reads.fastq`} />

This matches only the exact sequence `GACGAGATTA`, not the alternate allele `GACGTGATTA`.

For simple DNA patterns like `AGATCGGAAGAGC` that contain no regex metacharacters, `-F` and the default regex mode produce the same result, but `-F` can be faster on large files.

## The Tn5 Mosaic End

Let's find reads containing the Tn5 transposase mosaic end sequence, commonly seen in ATAC-seq and tagmentation-based library preps:

<Execute command={`fqgrep -c 'AGATGTGTATAAGAGACAG' reads.fastq`} />

<Execute command={`fqgrep 'AGATGTGTATAAGAGACAG' reads.fastq | head -8`} />

> **Tip**: The `-c` flag is your friend for quickly assessing contamination levels or checking how many reads contain a sequence of interest.
