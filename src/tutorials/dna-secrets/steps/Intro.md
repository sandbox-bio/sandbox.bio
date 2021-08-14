<script>
import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

<Alert>
	**Note**: This problem set assumes you already completed the tutorials on <Link href="/tutorials?id=bowtie2-intro">bowtie2</Link>, <Link href="/tutorials?id=samtools-intro">samtools</Link>, and <Link href="/tutorials?id=bedtools-intro">bedtools</Link>, but googling works too 🙂
</Alert>

In this problem set, you're given DNA sequencing reads that encode a secret message. To decode this message, you'll have to align those reads to a reference genome and call variants. The SNPs you'll identify via variant calling, when sorted by genomic coordinate, make up a string of DNA letters (`A`, `C`, `G`, `T`) that we'll decode to reveal the secret message.
