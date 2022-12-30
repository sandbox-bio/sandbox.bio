<script>
import Execute from "components/Execute.svelte";
import Exercise from "components/Exercise.svelte";

const criteria = [{
	name: "File <code>chromosomes.sorted.csv</code> contains a sorted version of file <code>chromosomes.csv</code>",
	checks: [{
		type: "file",
		path: "chromosomes.sorted.csv",
		action: "contents",
		commandExpected: `sort -V -k1,1 chromosomes.csv`
	}]
}];

const hints = [
];
</script>

The file `chromosomes.csv` is a CSV file (comma-separated values) that contains an **unsorted** list of chromosome names in the first column, and their size in the second column:

<Execute command={"cat chromosomes.csv"} />

We can use the `sort` command to sort this file by the first column, where `-n` signifies to sort numerically, and `-k1,1` means to only sort on the first column

<Execute command={"sort -n -k1,1 chromosomes.csv"} />

While the file is more sorted than before, it doesn't look quite right: `chr10` appears right after `chr1` but before `chr2`.

**Your goal**: Find the flag that allows you to sort the first column alphanumerically, i.e. `chr1`, `chr2`, `chr3`, ..., `chr10`. One done, save the results to the file `chromosomes.sorted.csv`

<Exercise {criteria} {hints} />
