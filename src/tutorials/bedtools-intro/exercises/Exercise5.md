<script>
// Solution:
//    bedtools intersect -a exons.bed -b <(grep Enhancer hesc.chromHmm.bed) -wa -wb -f 1.0 | wc -l > count

import Alert from "components/Alert.svelte";
import Exercise from "components/Exercise.svelte";

let criteria = [
{
	name: "File <code>count</code> contains the number of exons completely overlapped by an enhancer",
	checks: [{
		type: "file",
		path: "count",
		action: "contents",
		command: "bedtools intersect -a exons.bed -b <(grep Enhancer hesc.chromHmm.bed) -wa -wb -f 1.0 | wc -l"
	}]
}
];
</script>

Are there any exons that are completely overlapped by an enhancer? If so, how many?

<Alert>
	**Hint**: to find enhancers, you can use `grep` on the `hesc.chromHmm.bed` file using the `Enhancer` pattern.
</Alert>

<Exercise {criteria} />
