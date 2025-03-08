<script>
// Solution:
//    bedtools jaccard -a cpg.bed -b <(grep Enhancer hesc.chromHmm.bed) > jaccard.enhancers.txt; bedtools jaccard -a cpg.bed -b <(grep Promoter hesc.chromHmm.bed) > jaccard.promoters.txt

import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>jaccard.enhancers.txt</code> contains Jaccard stats between CpG and enhancers",
	checks: [{
		type: "file",
		path: "jaccard.enhancers.txt",
		action: "contents",
		commandExpected: "bedtools jaccard -a cpg.bed -b <(grep Enhancer hesc.chromHmm.bed)"
	}]
},
{
	name: "File <code>jaccard.promoters.txt</code> contains Jaccard stats between CpG and promoters",
	checks: [{
		type: "file",
		path: "jaccard.promoters.txt",
		action: "contents",
		commandExpected: "bedtools jaccard -a cpg.bed -b <(grep Promoter hesc.chromHmm.bed)"
	}]
},
];
</script>

What is the Jaccard statistic between CpG and hESC enhancers? Compare that to the Jaccard statistic between CpG and hESC promoters. Does the result make sense?

> **Hint**: to find enhancers, you can use `grep` on the `hesc.chromHmm.bed` file using the `Enhancer` pattern.

<Exercise {criteria} />
