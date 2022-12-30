<script>
import Execute from "components/Execute.svelte";
import Exercise from "components/Exercise.svelte";

const criteria = [{
	name: "File <code>exons.bed</code> no longer causes <code>bedtools</code> error <code>unable to determine types</code>",
	checks: [{
		type: "file",
		path: "exons.bed",
		action: "contents",
		commandExpected: `sed 's/ /\t/g' exons.bed`
	}]
}];
</script>

You just received a file from your collaborator: `exons.bed`, which contains a list of exonic regions. As part of your analysis, you would like to merge all overlapping regions into contiguous intervals.

This is a perfect job for `bedtools merge`, **but** you keep getting an error when you run the following command:

<Execute command={"bedtools merge -i exons.bed"} />

Your goal for this puzzle is to fix what is wrong with exons.bed such that you are able to run the command above without error.

<Exercise {criteria} />
