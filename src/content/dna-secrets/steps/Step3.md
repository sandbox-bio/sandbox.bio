<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.bam aligned.sam;  bcftools mpileup -f $REF_FASTA aligned.bam | bcftools call -m -v -Ob -o variants.bcf -; bcftools index variants.bcf
*/

import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>variants.bcf</code> exists",
	checks: [{
		type: "file",
		path: "variants.bcf",
		action: "exists",
	}]
},
{
	name: "File <code>variants.bcf</code> contains the variants",
	checks: [{
		type: "file",
		path: "variants.bcf",
		action: "contents",
		commandExpected: "bcftools mpileup -f $REF_FASTA aligned.bam | bcftools call -m -v -Ob - -o /shared/tmp/__dnasecret.bcf; bcftools view --no-header /shared/tmp/__dnasecret.bcf",
		commandObserved: "bcftools view --no-header variants.bcf"
	}]
},
{
	name: "File <code>variants.bcf.csi</code> is the index file of <code>variants.bcf</code> obtained with <code>bcftools index</code>",
	checks: [{
		type: "file",
		path: "variants.bcf.csi",
		action: "exists"
	}]
},
];
</script>

Now that we have reads aligned to the reference genome, let's call variants using `bcftools`. Output the variants to the file `variants.bcf`.

<Alert>
	**Hint**: Check out the <Link href="/tutorials/bowtie2-intro/6">bcftools section</Link> of the bowtie2 tutorial for an example of how to run `bcftools`.
</Alert>

<Exercise {criteria} />
