<script>
import { Icon } from "sveltestrap";
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

In this tutorial we'll analyze FASTA sequences of microRNA hairpins from the <Link href="https://www.mirbase.org/download/">miRNA database</Link>, and FASTQ sequencing reads from <Link href="https://42basepairs.com/browse/r2/genomics-data?file=reads_NA12878_R1.fastq.gz">42basepairs</Link>.

The data is preloaded here as `hairpins.fa` and `NA12878.fastq`: <Execute inline command='ls' />

Let's use SeqKit to calculate summary statistics for these two files:

<Execute command={`seqkit stats hairpins.fa NA12878.fastq`} />

<Alert color="primary">
    <Icon name="question-circle-fill" /> How did SeqKit know that `hairpins.fa` contains RNA sequences?
</Alert>

To avoid manually writing each file name, we can use a wildcard (`*`) to analyze all `.fa` and `.fastq` files in the current folder:

<Execute command={`seqkit stats *.{fa,fastq}`} />

<Alert>
    <Icon name="lightbulb-fill" /> Wildcards are a broadly useful feature on the command line.

    For example, `ls *.fa` lists all files ending in the extension `.fa` whereas <code>ls *.&lbrace;a,b,c}</code> lists files ending in either extension `.a`, `.b`, or `.c`.

</Alert>

<hr />

SeqKit can also calculate additional stats such as GC content, and the fraction of FASTQ reads with a mapping quality of 30. To enable those stats, use the `--all` flag:

<Execute command={`seqkit stats *.{fa,fastq} --all`} />

<Alert>
    <Icon name="lightbulb-fill" /> The columns are wrapped to the next line, which makes the output difficult to read.

    Try appending `| less -S` to the end of the command above, and use the right/left arrows to navigate columns without wrapping.

    Use `q` to exit `less`.

</Alert>
