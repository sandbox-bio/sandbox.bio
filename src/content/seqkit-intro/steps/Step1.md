<script>
import { Icon } from "sveltestrap";
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

In this tutorial we'll analyze FASTA sequences of microRNA hairpins from the <Link href="https://www.mirbase.org/download/">miRNA database</Link>, and FASTQ sequencing reads from <Link href="https://42basepairs.com/browse/r2/genomics-data?file=reads_NA12878_R1.fastq.gz">42basepairs</Link>.

The data is preloaded here as `hairpins.fa` and `NA12878.fastq`: <Execute inline command='ls' />

Let's use SeqKit to calculate summary statistics for these hairpin sequences:

<Execute command={`seqkit stats *`} />

<Alert color="primary">
    <Icon name="question-circle-fill" /> How did SeqKit know that `hairpins.fa` contains RNA sequences?
</Alert>
