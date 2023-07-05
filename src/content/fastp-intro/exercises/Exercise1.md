<script>
// Solution:
//    fastp --in1 HG004_R1.fastq.gz --in2 HG004_R2.fastq.gz --reads_to_process 10000

import Link from "$components/Link.svelte";
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "You ran <code>fastp</code> on only the first 10,000 reads",
	checks: [{
		type: "file",
		path: "fastp.json",
		action: "contents",
		commandExpected: "fastp --in1 HG004_R1.fastq.gz --in2 HG004_R2.fastq.gz --reads_to_process 10000 --html /shared/tmp/__fastp.html --json /shared/tmp/__fastp.json; jq '.summary' /shared/tmp/__fastp.json",
		commandObserved: "jq '.summary' fastp.json"
	}]
}];
</script>

Sometimes, a sequencing run generates a lot of data. Say we're only interested in getting a quick preview of data quality without having to analyze our entire dataset.

Use the <Link href="https://github.com/OpenGene/fastp/tree/v0.20.1#all-options">fastp README</Link> to find a parameter that allows you to **only process the first 10,000 reads** of the FASTQ files `HG004_R1.fastq.gz` and `HG004_R2.fastq.gz`.

<Exercise {criteria} />
