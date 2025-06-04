<script>
import Execute from "$components/Execute.svelte";
</script>

We now have SNPs in `variants.bcf` and `variants2.bcf`.

We want those SNPs to end up in the right genomic order, so let's use the `bcftools merge` command:

<Execute command="bcftools merge variants.bcf variants2.bcf > combined.vcf" />

Inspect the values in the `POS` column to make sure they are in the right order:

<Execute command="bcftools view combined.vcf" />
