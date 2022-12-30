<script>
import Execute from "components/Execute.svelte";
import Exercise from "components/Exercise.svelte";

const criteria = [{
	name: "File <code>chromosomes.sorted.csv</code> contains a sorted version of file <code>chromosomes.csv</code>",
	checks: [{
		type: "file",
		path: "chromosomes.sorted.csv",
		action: "contents",
		commandExpected: `sort -V chromosomes.csv`
	}]
}];

const hints = [
	"When googling for command line utilities like <code>sort</code>, prefix queries with <code>unix</code> or <code>linux</code> so it's clear which context you're interested in.",
	"Try googling <code>unix sort alphanumerical</code>."
];
</script>

The file `chromosomes.csv` is a CSV file (comma-separated values) that contains an **unsorted** list of chromosome names in the first column, and their size in the second column:

<Execute command={"cat chromosomes.csv"} />

We can use the `sort` command to sort this file, where `-n` signifies to sort numerically *(usually you'll want to specify which columns to sort on using `-k`, but for this exercise, it's ok to sort on all columns)*:

<Execute command={"sort -n chromosomes.csv"} />

While the file is more sorted than before, it doesn't look quite right: `chr10` appears right after `chr1` but before `chr2`.

**Your goal**: Find the flag that allows you to sort the file alphanumerically, i.e. `chr1`, `chr2`, `chr3`, ..., `chr10`. One done, save the results to the file `chromosomes.sorted.csv`.

<Exercise {criteria} {hints} />
