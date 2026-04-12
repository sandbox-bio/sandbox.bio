<script>
import Execute from "$components/Execute.svelte";
</script>

In the previous step, we manually translated IUPAC codes like `N` into regex character classes like `[ACGT]`. fqgrep's `--iupac` flag automates this — it requires `-F` (fixed-string mode) and converts IUPAC codes in your pattern before matching.

## The --iupac flag

The `--iupac` flag has four modes that control how IUPAC ambiguity codes are handled during fixed-string matching:

| Mode       | Behavior                           |
| ---------- | ---------------------------------- |
| `never`    | Default — patterns are left as-is  |
| `expand`   | Expand into multiple fixed strings |
| `regex`    | Convert to regex character classes |
| `bit-mask` | Use 4-bit bitmask encoding         |

## Expand mode

With `--iupac expand`, each IUPAC code is expanded into all possible fixed strings. For example, `K` (G or T) in `GATK` produces two patterns: `GATG` and `GATT`.

Let's search for our degenerate barcode `ACGTNNNNTGCA` using expand mode:

<Execute command={`fqgrep -F --iupac expand -c 'ACGTNNNNTGCA' reads.fastq`} />

Each `N` is expanded to A, C, G, T — producing all possible concrete sequences. This can get large: 4 N's produce 4^4 = 256 patterns!

## Regex mode

With `--iupac regex`, IUPAC codes are converted to regex character classes. This is often more efficient than `expand` for patterns with many degenerate positions:

<Execute command={`fqgrep -F --iupac regex -c 'ACGTNNNNTGCA' reads.fastq`} />

Behind the scenes, `ACGTNNNNTGCA` becomes `ACGT[ACGT][ACGT][ACGT][ACGT]TGCA` — exactly what we wrote by hand in the previous step!

## Bit-mask mode

The `bit-mask` mode uses a compact 4-bit encoding where each IUPAC code is the bitwise OR of its constituent bases (A=1, C=2, G=4, T=8):

<Execute command={`fqgrep -F --iupac bit-mask -c 'ACGTNNNNTGCA' reads.fastq`} />

This mode can be efficient for patterns with many degenerate positions.

## Common IUPAC codes

Here's a quick reference for the most common IUPAC ambiguity codes:

| Code | Bases      | Meaning    |
| ---- | ---------- | ---------- |
| N    | A, C, G, T | Any base   |
| R    | A, G       | Purine     |
| Y    | C, T       | Pyrimidine |
| K    | G, T       | Keto       |
| M    | A, C       | Amino      |
| S    | C, G       | Strong     |
| W    | A, T       | Weak       |

## When to use each mode

- **`expand`**: Best for patterns with few degenerate positions (1-3 N's). Simple and fast.
- **`regex`**: Best for patterns with many degenerate positions. Avoids combinatorial explosion.
- **`bit-mask`**: Efficient alternative to regex for heavily degenerate patterns.

> **Key takeaway**: The `--iupac` flag lets you search for degenerate sequences naturally using standard IUPAC notation, without manually translating to regex. Use `-F --iupac regex` for most cases.
