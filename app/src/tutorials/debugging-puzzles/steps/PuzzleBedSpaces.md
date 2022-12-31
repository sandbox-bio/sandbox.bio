<script>
import Execute from "components/Execute.svelte";
import Exercise from "components/Exercise.svelte";

const criteria = [{
	name: "File <code>exons.merged.bed</code> no longer causes <code>bedtools merge</code> to output an error",
	checks: [{
		type: "file",
		path: "exons.merged.bed",
		action: "contents",
		commandExpected: `sed 's/ /\t/g' exons.bed`
	}]
}];

const hints = [
	`As suggested by the error message, run <code>cat -t exons.bed</code>. Also try <code>head exons.bed</code>. Do any lines stand out from the others?`,
	`In the output of <code>cat -t exons.bed</code>, the first line uses spaces as the column delimiter instead of tabs.`,
	`You can use a <code>sed</code> command to replace spaces with tabs (<code>\\t</code>).`,
	`Don't forget to specify that you want the <code>sed</code> replace logic to be global!`
];
</script>

You just received a file from your collaborator: `exons.bed`, which contains a list of exonic regions. As part of your analysis, you would like to merge all overlapping regions into contiguous intervals.

This is a perfect job for `bedtools merge`, **but** you keep getting an error when you run the following command:

<Execute command={"bedtools merge -i exons.bed"} />

**Your Goal**: Create a file `exons.merged.bed` from the regions in `exons.bed` so that the following query does not give an error:

<Execute command={"bedtools merge -i exons.merged.bed"} />

<Exercise {criteria} {hints} />
