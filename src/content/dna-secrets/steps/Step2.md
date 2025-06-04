<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.bam aligned.sam
*/

import Link from "$components/Link.svelte";
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>aligned.sam</code> contains reads mapped to the genome using <code>bowtie2</code>",
	checks: [{
		type: "file",
		path: "aligned.sam",
		action: "contents",
		commandExpected: `echo -e "303 aligned.sam"`,
		commandObserved: "wc -l aligned.sam"
	}]
},
{
	name: "File <code>aligned.bam</code> is a sorted BAM file version of <code>aligned.sam</code>",
	checks: [{
		type: "file",
		path: "aligned.bam",
		action: "contents",
		commandExpected: "samtools sort -o /tmp/__dnasecret.bam aligned.sam; samtools view /tmp/__dnasecret.bam",
		commandObserved: "samtools view aligned.bam",
	}]
}];
</script>

First, use `bowtie2` to align the sequencing reads in `reads.fq` to the reference genome using the index located at `$REF`; the reads are single-ended.

> **Hint**: Check out the <Link href="/tutorials/bowtie2-intro/2">bowtie2 tutorial</Link> for an example of how to use `bowtie2`.

Output the results to the file `aligned.sam`, then sort the SAM file to output `aligned.bam`. Complete the following exercises before moving on to the next step:

<Exercise {criteria} />
