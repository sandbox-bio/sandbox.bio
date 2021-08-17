<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.sorted.bam aligned.sam; bcftools mpileup -f $REF_FASTA aligned.sorted.bam | bcftools call -m -v -Ov -o variants.vcf -; bcftools query --format "%POS\t%ALT\n" variants.vcf > secret.tsv
*/

import Execute from "./components/Execute.svelte";
</script>

Oups. I forgot to mention that we have a second set of DNA sequencing data, stored in the file `morereads.fq` 😬

The secret you obtained in `secret.tsv` is only half the puzzle, so we'll need to call variants on the file `morereads.fq`, and merge our results in the correct genomic order.

To save us some time, here's a very long pipeline you can use to perform read alignment and variant calling:

<Execute command='bowtie2 \ -x $REF \ -U morereads.fq \ -S aligned2.sam; \ samtools sort -o aligned2.sorted.bam aligned2.sam; \ bcftools mpileup -f $REF_FASTA aligned2.sorted.bam | \ bcftools call -m -v -Ov -o variants2.vcf -; \ bcftools query --format "%POS\t%ALT\\n" variants2.vcf > secret2.tsv' />

In the next step, we'll see how we can merge those `secrets.tsv` and `secrets2.tsv` in the correct genomic order.
