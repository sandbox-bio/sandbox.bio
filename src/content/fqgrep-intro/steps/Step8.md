<script>
import Execute from "$components/Execute.svelte";
</script>

fqgrep supports full regular expressions, making it powerful for searching degenerate sequences, variable regions, or complex patterns.

## Basic regex patterns

Let's find reads containing a barcode with a variable middle section. The pattern `ACGT....TGCA` matches `ACGT`, followed by any 4 bases, followed by `TGCA`:

<Execute command={`fqgrep -c 'ACGT....TGCA' reads.fastq`} />

<Execute command={`fqgrep --color always 'ACGT....TGCA' reads.fastq | head -16`} />

The `.` matches any single character (nucleotide).

## Character classes

Use `[ACGT]` to match any of those specific characters:

<Execute command={`fqgrep --color always 'ACGT[ACGT][ACGT][ACGT][ACGT]TGCA' reads.fastq | head -8`} />

This is equivalent to `ACGT....TGCA` but more explicit about matching only nucleotides.

## IUPAC ambiguity codes

In molecular biology, IUPAC codes represent degenerate bases:

- `N` = any base (A, C, G, or T)
- `R` = purine (A or G)
- `Y` = pyrimidine (C or T)
- `W` = weak (A or T)
- `S` = strong (C or G)

You can express these as regex character classes:

| IUPAC | Regex    | Meaning    |
| ----- | -------- | ---------- |
| N     | `[ACGT]` | Any base   |
| R     | `[AG]`   | Purine     |
| Y     | `[CT]`   | Pyrimidine |
| W     | `[AT]`   | Weak       |
| S     | `[CG]`   | Strong     |

## Example: Degenerate primer

Imagine searching for a primer with sequence `ACGTNNNNTGCA` where `N` means any base:

<Execute command={`fqgrep --color always 'ACGT[ACGT][ACGT][ACGT][ACGT]TGCA' reads.fastq | head -12`} />

Translating IUPAC codes to regex by hand works, but fqgrep also has a dedicated `--iupac` flag to do this automatically — we'll cover that in the next step.

## Quantifiers

Use `{n}` to specify exact repetitions:

- `A{5}` matches exactly 5 A's
- `[ACGT]{8}` matches any 8 nucleotides

Let's find reads with poly-A stretches (10 or more A's):

<Execute command={`fqgrep -c 'A{10}' reads.fastq`} />

<Execute command={`fqgrep --color always 'A{10}' reads.fastq | head -8`} />

## Alternation

Use `|` to match alternative patterns:

<Execute command={`fqgrep -c 'GACGAGATTA|GACGTGATTA' reads.fastq`} />

This finds reads with either the reference (`GACGAGATTA`) or alternate (`GACGTGATTA`) allele — useful for variant detection!

> **Tip**: When using regex special characters in your shell, always quote your patterns to prevent shell interpretation.
