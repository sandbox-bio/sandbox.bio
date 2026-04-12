<script>
import Execute from "$components/Execute.svelte";
</script>

fqgrep isn't limited to DNA — the `--protein` flag lets you search FASTQ files containing protein (amino acid) sequences.

## The --protein flag

When `--protein` is specified, fqgrep treats input sequences as amino acids instead of DNA. This flag cannot be combined with DNA-specific options like `--reverse-complement` or `--iupac` — fqgrep will report an error if you try.

Let's look at our small protein FASTQ file:

<Execute command={`head -12 proteins.fastq`} />

Notice the sequences contain amino acid characters (M, K, Y, L, P, etc.) rather than just A, C, G, T.

## Searching for protein motifs

Search for a signal peptide motif:

<Execute command={`fqgrep --protein -c 'MKYLLPTAAAGLLLLAAQ' proteins.fastq`} />

<Execute command={`fqgrep --protein --color always 'MKYLLPTAAAGLLLLAAQ' proteins.fastq | head -8`} />

## Finding His tags

Polyhistidine tags (His tags) are commonly used in protein purification. Let's find reads containing a 6xHis tag:

<Execute command={`fqgrep --protein -c 'HHHHHH' proteins.fastq`} />

## Regex with protein sequences

Regular expressions work with protein sequences too:

<Execute command={`fqgrep --protein --color always 'H{6}' proteins.fastq | head -8`} />

## Limitations

The `--protein` flag cannot be used with:

- `--reverse-complement` (DNA-specific concept)
- `--iupac` (DNA ambiguity codes, not amino acid codes)

> **Tip**: The `--protein` flag is useful for searching FASTQ output from protein sequencing platforms or translated sequence databases stored in FASTQ format.
