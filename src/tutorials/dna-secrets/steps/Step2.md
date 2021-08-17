<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.sorted.bam aligned.sam
*/

import Exercise from "./components/Exercise.svelte";

let criteria = [
{
	name: "File <code>aligned.sam</code> contains reads mapped to the genome using <code>bowtie2</code>",
	checks: [{
		type: "file",
		path: "aligned.sam",
		action: "contents",
		commandExpected: "bowtie2 -x $REF -U /shared/data/reads.fq | samtools view",
		commandObserved: "samtools view aligned.sam"
	}]
},
{
	name: "File <code>aligned.sorted.bam</code> is a sorted BAM file version of <code>aligned.sam</code>",
	checks: [{
		type: "file",
		path: "aligned.sorted.bam",
		action: "contents",
		commandExpected: "samtools sort -o /shared/tmp/__dnasecret.bam aligned.sam; samtools view /shared/tmp/__dnasecret.bam",
		commandObserved: "samtools view aligned.sorted.bam",
	}]
}];
</script>

First, use `bowtie2` to align the sequencing reads in `reads.fq` to the reference genome using the index located at `$REF`; the reads are single-ended. Output the resulting SAM file to the file `aligned.sam`.

Then, sort the SAM file and index it. Complete the following exercises before moving on to the next step:

<Exercise {criteria} />
