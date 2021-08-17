<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.sorted.bam aligned.sam;  bcftools mpileup -f $REF_FASTA aligned.sorted.bam | bcftools call -m -v -Ob -o variants.bcf -; bcftools index variants.bcf

	bowtie2 -x $REF -U morereads.fq -S aligned2.sam; samtools sort -o aligned2.sorted.bam aligned2.sam;  bcftools mpileup -f $REF_FASTA aligned2.sorted.bam | bcftools call -m -v -Ob -o variants2.bcf -; bcftools index variants2.bcf
*/

import Execute from "./components/Execute.svelte";
</script>

At this point, we have SNPs in `variants.bcf` and `variants2.bcf` across the genome. We need to merge those two BCF files so they end up in the right genomic order.

To do so, simply use the `bcftools merge` command:

<Execute command="bcftools merge variants.bcf variants2.bcf" />

Inspect the values in the `POS` column to make sure they are in the right order.

To see a simpler view with just the position, reference allele, and alternate allele, use `bcftools query`:

<Execute command='bcftools merge variants.bcf variants2.bcf | \ bcftools query -f "%POS\t%REF\t%ALT\n" -' />

Let's now see if we can finally decode this DNA secret
