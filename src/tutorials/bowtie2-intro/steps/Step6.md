<script>
import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

To generate variant calls in the VCF format, run:

<Execute command="bcftools mpileup -f $REF_FA eg2.sorted.bam | \ bcftools view -Ov -o eg2.raw.bcf -" />

<Alert>
	**Note**: We preloaded your sandbox with the variable `$REF_FA`. You can use the <Execute command="env" inline={true} /> command to list variables and <Execute command="echo $REF_FA" inline={true} /> to view their value.
</Alert>

Then to view the variants, run:

<Execute command="bcftools view eg2.raw.bcf | head" />

See the official bcftools <Link href="http://samtools.github.io/bcftools/howtos/variant-calling.html">guide to variant calling</Link> for more details and variations on this process.
