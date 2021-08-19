<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.sorted.bam aligned.sam;  bcftools mpileup -f $REF_FASTA aligned.sorted.bam | bcftools call -m -v -Ob -o variants.bcf -; bcftools index variants.bcf
*/

import Execute from "./components/Execute.svelte";
</script>

Now that we're done with variant calling, let's see what these variants look like:

<Execute command="bcftools view variants.bcf" />

The lines that start with `#` are header lines that contain metadata about the VCF file. To remove those lines:

<Execute command="bcftools view --no-header variants.bcf" />

There's a lot of information in VCF files but for our purposes we don't need most of it, so let's use `bcftools query` to subset our output. For example, let's only show the genomic position, the reference allele and the allele that was called at that position:

<Execute command='bcftools query --format "%POS\t%REF\t%ALT\n" variants.bcf' />

The last column (`%ALT`) contains the SNPs we called, and correspond to the secret message that is encoded into DNA.
