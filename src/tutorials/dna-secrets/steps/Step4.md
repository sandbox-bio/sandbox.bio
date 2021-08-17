<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.sorted.bam aligned.sam; bcftools mpileup -f $REF_FASTA aligned.sorted.bam | bcftools call -m -v -Ov -o variants.vcf -; bcftools query --format "%POS\t%ALT\n" variants.vcf > secret.tsv
*/

import Execute from "./components/Execute.svelte";
import Exercise from "./components/Exercise.svelte";

let criteria = [
{
	name: "File <code>secret.tsv</code> is a tab-separated file with 2 columns, the genomic coordinate and the SNP at that position",
	checks: [{
		type: "file",
		path: "secret.tsv",
		action: "contents",
		command: 'bcftools query --format "%POS\t%ALT\n" variants.vcf'
	}]
}];
</script>

Now that we're done with variant calling, let's see what those variants look like:

<Execute command="bcftools view variants.vcf" />

The lines that start with `#` are header lines that contain metadata about the VCF file. To remove those lines:

<Execute command="bcftools view --no-header variants.vcf" />

There's a lot of information in VCF files but for our purposes we don't need most of it, so let's use `bcftools query` to subset our output. For example, let's only show the genomic position, the reference allele and the allele that was called at that position:

<Execute command='bcftools query --format "%POS\t%REF\t%ALT\n" variants.vcf' />

The last column (`%ALT`) are the SNPs we called, and correspond to the secret message that is encoded into DNA.

Let's store that information for later: Create a file `secret.tsv`, a tab-separated file with two columns: `%POS` and `%ALT`:

<Exercise {criteria} />
