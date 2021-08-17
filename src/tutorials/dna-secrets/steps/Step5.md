<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.sorted.bam aligned.sam;  bcftools mpileup -f $REF_FASTA aligned.sorted.bam | bcftools call -m -v -Ob -o variants.bcf -; bcftools index variants.bcf

	bowtie2 -x $REF -U morereads.fq -S aligned2.sam; samtools sort -o aligned2.sorted.bam aligned2.sam;  bcftools mpileup -f $REF_FASTA aligned2.sorted.bam | bcftools call -m -v -Ob -o variants2.bcf -; bcftools index variants2.bcf
*/

import Execute from "./components/Execute.svelte";
</script>

Oups. I forgot to mention that we have a second set of DNA sequencing data, stored in the file `morereads.fq` 😬

The SNPs you obtained in `variants.bcf` are only half the puzzle; we'll need to call variants on the file `morereads.fq` and merge the VCFs in the correct genomic order so that we can decode the secret message.

To save us some time, here's a pipeline you can use to perform read alignment and variant calling on that second dataset:

<Execute command='bowtie2 -x $REF \ -U morereads.fq \ -S aligned2.sam; \ samtools sort -o aligned2.sorted.bam aligned2.sam; \ bcftools mpileup -f $REF_FASTA aligned2.sorted.bam | \ bcftools call -m -v -Ob -o variants2.bcf -; bcftools index variants2.bcf' />

In the next step, we'll see how to merge `variants.bcf` and `variants2.bcf` in the correct genomic order.
