<script>
// Solution:
//    grep Enhancer hesc.chromHmm.bed > enhancers.bed
//    grep Promoter hesc.chromHmm.bed > promoters.bed
//    bedtools jaccard -a cpg.bed -b enhancers.bed > jaccard.enhancers.txt
//    bedtools jaccard -a cpg.bed -b promoters.bed > jaccard.promoters.txt

// import { CoreUtils } from "terminal/coreutils";
import Exercise from "components/Exercise.svelte";
import Alert from "components/Alert.svelte";

let criteria = [
{
	name: "File <code>enhancers.bed</code> contains a list of enhancers",
	checks: [{
		type: "file",
		path: "enhancers.bed",
		action: "contents",
		fn: async () => await CoreUtils.grep(["Enhancer", "hesc.chromHmm.bed"]),
		output: "/shared/tmp/exercise4-enhancers.bed"
	}]
},
{
	name: "File <code>promoters.bed</code> contains a list of promoters",
	checks: [{
		type: "file",
		path: "promoters.bed",
		action: "contents",
		fn: async () => await CoreUtils.grep(["Promoter", "hesc.chromHmm.bed"]),
		output: "/shared/tmp/exercise4-promoters.bed"
	}]
},
{
	name: "File <code>jaccard.enhancers.txt</code> contains Jaccard stats between CpG and enhancers",
	checks: [{
		type: "file",
		path: "jaccard.enhancers.txt",
		action: "contents",
		command: "bedtools jaccard -a cpg.bed -b /shared/tmp/exercise4-enhancers.bed"
	}]
},
{
	name: "File <code>jaccard.promoters.txt</code> contains Jaccard stats between CpG and promoters",
	checks: [{
		type: "file",
		path: "jaccard.promoters.txt",
		action: "contents",
		command: "bedtools jaccard -a cpg.bed -b /shared/tmp/exercise4-promoters.bed"
	}]
},
];
</script>

What is the Jaccard statistic between CpG and hESC enhancers? Compare that to the Jaccard statistic between CpG and hESC promoters. Does the result make sense?

<Alert>
	**Hint**: to find enhancers, you can use `grep` on the `hesc.chromHmm.bed` file using the `Enhancer` pattern.
</Alert>

<Exercise {criteria} />
