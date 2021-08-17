<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.sorted.bam aligned.sam; bcftools mpileup -f $REF_FASTA aligned.sorted.bam | bcftools call -m -v -Ov -o variants.vcf -
*/

import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";
import Exercise from "./components/Exercise.svelte";

let criteria = [
{
	name: "File <code>variants.vcf</code> contains variants called using <code>bcftools</code> (make sure to <strong>output a VCF file</strong>, not a BCF)",
	checks: [{
		type: "file",
		path: "variants.vcf",
		action: "contents",
		command: "bcftools mpileup -f $REF_FASTA aligned.sorted.bam | bcftools view -Ov -",
		filter: d => d.split("\n").filter(l => !l.startsWith("#")).join("\n")
	}]
}];
</script>

Now that we have reads aligned to the reference genome, let's call variants using `bcftools`. Output the variants to the file `variants.vcf`.

<Alert>
	**Hint**: Check out the <Link href="/tutorials?id=bowtie2-intro&step=6">bcftools section</Link> of the bowtie2 tutorial for an example of how to run `bcftools`.
</Alert>

<Exercise {criteria} />
