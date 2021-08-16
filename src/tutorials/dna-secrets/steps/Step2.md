<script>
// Solution:
//    bowtie2 -x $REF -U reads.fq -S aligned.sam
//    samtools sort -o aligned.sorted.bam aligned.sam


import Exercise from "./components/Exercise.svelte";

let criteria = [
{
	name: "File <code>aligned.sam</code> contains reads mapped to the genome",
	checks: [{
		type: "file",
		path: "aligned.sam",
		action: "contents",
		contents: "yes"
	}]
},
{
	name: "File <code>aligned.sorted.bam</code> is a sorted BAM file version of <code>aligned.sam</code>",
	checks: [{
		type: "file",
		path: "aligned.sorted.bam",
		action: "contents",
		command: "samtools sort -o /tmp/__dnasecret.bam aligned.sam; cat /tmp/__dnasecret.bam"
	}]
},
];
</script>

First, use `bowtie2` to align the sequencing reads in `reads.fq` to the reference genome using the index located at `$REF`. Then, sort the SAM file and index it.

Assume the reads are single-ended. Output the resulting SAM file to the file `aligned.sam`.

<Exercise {criteria} />
