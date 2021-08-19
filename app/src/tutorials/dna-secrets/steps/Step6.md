<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.bam aligned.sam;  bcftools mpileup -f $REF_FASTA aligned.bam | bcftools call -m -v -Ob -o variants.bcf -; bcftools index variants.bcf

	bowtie2 -x $REF -U morereads.fq -S aligned2.sam; samtools sort -o aligned2.bam aligned2.sam;  bcftools mpileup -f $REF_FASTA aligned2.bam | bcftools call -m -v -Ob -o variants2.bcf -; bcftools index variants2.bcf

	bcftools merge variants.bcf variants2.bcf > combined.vcf
*/

import Execute from "./components/Execute.svelte";
</script>

We now have SNPs in `variants.bcf` and `variants2.bcf`.

We want those SNPs to end up in the right genomic order, so let's use the `bcftools merge` command:

<Execute command="bcftools merge variants.bcf variants2.bcf > combined.vcf" />

Inspect the values in the `POS` column to make sure they are in the right order:

<Execute command="bcftools view combined.vcf" />
