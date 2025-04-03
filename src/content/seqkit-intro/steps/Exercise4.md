<script>
// Solution:
//    seqkit replace --pattern 'SRR000926.([0-9]+)' --replacement 'NA12878_read$1' NA12878.fastq > renamed.fastq

import Link from "$components/Link.svelte";
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>renamed.fastq</code> exists",
	checks: [{
		type: "file",
		path: "renamed.fastq",
		action: "exists"
	}]
},
{
	name: "File <code>renamed.fastq</code> contains sequences from <code>hairpins.fa</code> with renamed IDs",
	checks: [
		{
		type: "file",
		path: "renamed.fastq",
		action: "contents",
        commandObserved: "echo 1000",
        commandExpected: "grep NA12878 renamed.fastq | wc -l"  // much faster check than running seqkit seq
	},
		{
		type: "file",
		path: "renamed.fastq",
		action: "contents",
        commandObserved: "echo 0",
        commandExpected: "grep SRR000926 renamed.fastq | wc -l"  // much faster check than running seqkit seq
	}
	]
}];
</script>

This exercise is more challenging as it involves regular expressions, but is common enough that it's useful to work through it.

The FASTQ file `NA12878.fastq` has read IDs such as `@SRR000926.25` and your goal is to rename those to `@NA12878_read25` (you can keep the rest of the read name intact).

Use <Link href="https://bioinf.shenwei.me/seqkit/usage/#replace">seqkit replace</Link> to do so, and output the result to `renamed.fastq`.

<Exercise {criteria} />
