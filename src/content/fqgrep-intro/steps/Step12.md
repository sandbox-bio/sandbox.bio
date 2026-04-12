<script>
import Execute from "$components/Execute.svelte";
</script>

fqgrep is designed for speed. Let's explore what makes it fast and how to tune it.

## Why fqgrep is fast

Traditional approaches to searching FASTQ files have problems:

| Method         | Issue                                        |
| -------------- | -------------------------------------------- |
| `grep -B1 -A2` | Very slow, context lines break on edge cases |
| `awk` scripts  | Complex, slow, error-prone                   |
| Python/Perl    | Flexible but slow for large files            |

fqgrep is written in Rust and uses:

- **FASTQ-aware parsing** — no context-line hacks needed
- **Multi-threaded searching** — parallel pattern matching across reads
- **Efficient regex engine** — optimized Aho-Corasick for fixed strings, fast DFA for regex

In [benchmarks](https://github.com/Rbfinch/grepq) searching 30 patterns across an 874MB FASTQ, fqgrep finished in **0.34 seconds** vs. grep's **344 seconds** — roughly 1000x faster.

## Threads

The `--threads` option controls how many threads fqgrep uses for pattern matching. By default, fqgrep uses as many threads as there are CPUs on your machine.

In this sandbox environment, only one CPU is available:

<Execute command={`fqgrep --threads 1 -c 'AGATCGGAAGAGC' reads.fastq`} />

On a machine with multiple CPUs, increasing threads speeds up searches on large files (millions of reads). The threading model works as follows:

- One or more **reader** threads decompress and read input
- Multiple **matcher** threads search reads for patterns in parallel
- One **writer** thread outputs results

## Output ordering

By default, fqgrep preserves the input order of records in the output, even when using multiple threads. This means your output FASTQ will have reads in the same order as the input.

If order doesn't matter and you want to reduce memory usage slightly, use `--no-order`:

```bash
fqgrep --no-order 'AGATCGGAAGAGC' reads.fastq
```

## Compression support

fqgrep automatically handles gzip-compressed files based on the `.gz` extension:

```bash
fqgrep 'AGATCGGAAGAGC' reads.fastq.gz
```

You can also force decompression with `-Z` (or `--decompress`) for files without the `.gz` extension.

## Fixed strings for speed

When your pattern contains no regex metacharacters, using `-F` (or `--fixed-strings`) can be faster because fqgrep uses the Aho-Corasick algorithm for fixed-string matching instead of the regex engine:

<Execute command={`fqgrep -F -c 'AGATCGGAAGAGC' reads.fastq`} />

This difference is most noticeable on large files with many patterns (via `-f`).

## Exit codes

fqgrep follows `grep` conventions:

- `0` — One or more matches found
- `1` — No matches found
- `>1` — An error occurred

This is useful in scripts:

```bash
fqgrep -c 'AGATCGGAAGAGC' reads.fastq > /dev/null 2>&1 && echo "Adapter contamination detected!"
```

> **Summary**: For large-scale analysis, fqgrep's speed makes tasks that would take hours with grep finish in seconds. The familiar grep-like interface means you can drop it into existing workflows easily.
