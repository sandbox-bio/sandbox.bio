<script>
import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

To generate variant calls in the VCF format, run `bcftools mpileup`, followed by `bcftools call`:

<Execute command="bcftools mpileup -f $REF_FASTA eg2.sorted.bam | \ bcftools call -m -v -Ob -o eg2.bcf -" />

where: 

* `-Ob`: we want the output to be a BCF file (i.e. binary VCF). To output a text VCF, use `-Ov`.
* `-m`: use the default variant calling method
* `-v`: only output the variants (i.e. don't list the sites where we match the reference sequence)

<Alert>
	**Note**: We preloaded your sandbox with the variable `$REF_FASTA`. You can use the <Execute command="env" inline={true} /> command to list variables and `echo` to view their value, e.g. <Execute command="echo $REF_FASTA" inline={true} />.
</Alert>

Then to view the variants, run:

<Execute command="bcftools view eg2.bcf" />

See the official bcftools <Link href="http://samtools.github.io/bcftools/howtos/variant-calling.html">guide to variant calling</Link> for more details and variations on this process.
