<script>
import Execute from "$components/Execute.svelte";
</script>

When inspecting matches visually, it helps to see exactly where in the sequence the pattern matched.

## Highlighting matches with --color

Use `--color always` to highlight the matching portion of each sequence:

<Execute command={`fqgrep --color always 'AGATCGGAAGAGC' reads.fastq | head -20`} />

The matched sequence will be highlighted, making it easy to see where in the read the adapter appears — often near the 3' end for adapter contamination.

Let's try with the Tn5 mosaic end:

<Execute command={`fqgrep --color always 'AGATGTGTATAAGAGACAG' reads.fastq | head -12`} />

## Progress reporting

For large files, you might want to see progress as fqgrep works. Use `--progress`:

<Execute command={`fqgrep --progress -c 'AGATCGGAAGAGC' reads.fastq`} />

With our small test file, this finishes instantly. On real sequencing files with millions of reads, progress reporting helps you estimate how long the search will take.

## Combining options

You can combine color highlighting with other options. Let's view colored output for barcodes:

<Execute command={`fqgrep --color always 'ACGTACGT' reads.fastq | head -12`} />

> **Note**: Use `--color always` when you want highlighting regardless of whether output goes to a terminal. Use `--color auto` (the default) to only highlight when outputting to a terminal.
