<script>
import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";
</script>

<Alert>
	**Note**: This problem set assumes you already completed the tutorials on <Link href="/tutorials?id=bowtie2-intro">bowtie2</Link>, <Link href="/tutorials?id=samtools-intro">samtools</Link>, and <Link href="/tutorials?id=bedtools-intro">bedtools</Link>, but googling works too 🙂
</Alert>

In this problem set, you're given DNA sequencing reads that encode a secret message. To _decode_ this message, you'll have to align these reads to a reference genome and call variants. The SNPs you identify via variant calling, when sorted by genomic coordinate, will make up a string of DNA letters (`A`, `C`, `G`, `T`) that can be decoded to reveal the secret message.
