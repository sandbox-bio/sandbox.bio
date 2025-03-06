<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

In this tutorial we'll analyze FASTA sequences of microRNA hairpins from the <Link href="https://www.mirbase.org/download">miRNA database</Link>, and FASTQ sequencing reads from <Link href="https://42basepairs.com/browse/r2/genomics-data?file=reads_NA12878_R1.fastq.gz">42basepairs</Link>.

The data is preloaded here as `hairpins.fa` and `NA12878.fastq`: <Execute inline command='ls' />

Let's use SeqKit to calculate summary statistics for these two files:

<Execute command={`seqkit stats hairpins.fa NA12878.fastq`} />

> How did SeqKit know that `hairpins.fa` contains RNA sequences?

To avoid manually writing each file name, we can use a wildcard (`*`) to analyze all `.fa` and `.fastq` files in the current folder:

<Execute command={`seqkit stats *.{fa,fastq}`} />

> 💡 Wildcards are a broadly useful feature on the command line.
>
> For example, `ls *.fa` lists all files ending in the extension `.fa` whereas `ls *.{fa,fastq}` lists files ending in **either** `.fa` or `.fastq`.

<hr />

SeqKit can also calculate additional stats such as GC content, and the fraction of FASTQ reads with a mapping quality of 30. To enable those stats, use the flag `--all`:

<Execute command={`seqkit stats *.{fa,fastq} --all`} />

> 💡 The columns are wrapped to the next line, which makes the output difficult to read.
>
> Try appending `| less -S` to the end of the command above, and use the right/left arrows to navigate columns without wrapping.
>
> Use `q` to exit `less`.
