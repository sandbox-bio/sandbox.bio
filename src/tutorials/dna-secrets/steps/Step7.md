<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.sorted.bam aligned.sam;  bcftools mpileup -f $REF_FASTA aligned.sorted.bam | bcftools call -m -v -Ob -o variants.bcf -; bcftools index variants.bcf

	bowtie2 -x $REF -U morereads.fq -S aligned2.sam; samtools sort -o aligned2.sorted.bam aligned2.sam;  bcftools mpileup -f $REF_FASTA aligned2.sorted.bam | bcftools call -m -v -Ob -o variants2.bcf -; bcftools index variants2.bcf
	
	bcftools merge variants.bcf variants2.bcf | bcftools query -f "%ALT" -o secret -
*/

import Execute from "./components/Execute.svelte";
</script>
