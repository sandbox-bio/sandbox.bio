<script>
import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

<Link href="http://www.htslib.org/">SAMtools</Link> is a collection of tools for manipulating and analyzing SAM and BAM alignment files.

Use `samtools view` to convert the SAM file from 2 steps ago into a BAM file. BAM is the binary format corresponding to the SAM text format:

<Execute command="samtools view eg2.sam -o eg2.bam" />

Use `samtools sort` to convert the BAM file to a sorted BAM file:

<Execute command="samtools sort eg2.sam -o eg2.sorted.bam" />

We now have a sorted BAM file called `eg2.sorted.bam`. Sorted BAM is a useful format because the alignments are (a) compressed, which is convenient for long-term storage, and (b) sorted, which is convenient for variant discovery.

<Alert>
	**Note**: To learn more about `samtools`, check out our <Link href="/tutorials?id=samtools-intro">samtools tutorial</Link>.
</Alert>

