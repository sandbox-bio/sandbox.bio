<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
import Exercise from "$components/Exercise.svelte";

const criteria = [{
	name: "File <code>chromosomes.sorted.txt</code> contains a sorted version of file <code>chromosomes.txt</code>",
	checks: [{
		type: "file",
		path: "chromosomes.sorted.txt",
		action: "contents",
		commandExpected: `sort -V chromosomes.txt`
	}]
}];

const hints = [
	"Browse the <code>sort</code> command's help page using <code>sort --help</code>.",
	"In that help page, look for a way to do a <code>natural sort</code>.",
	"The command <code>sort --help | grep natural</code> should help you find the flag you're looking for."
];
</script>

The file `chromosomes.txt` is a text file that contains an **unsorted** list of chromosome names:

<Execute command={"cat chromosomes.txt"} />

We can use the `sort` command to sort this file, where `-n` sorts numerically:

<Execute command={"sort -n chromosomes.txt"} />

While the file is more sorted than before, it doesn't look quite right: `chr10` appears right after `chr1` but before `chr2`.

**Your goal**: Find the `sort` flag that allows you to sort the file using a **natural sort**, i.e. `chr1`, `chr2`, `chr3`, ..., `chr10`. Save the sorted file as `chromosomes.sorted.txt`.

<Alert>
	**Natural sort** is useful in genomics for sorting chromosome names, but it's also used to sort software versions. For example, with versions `v2` and `v11`, natural sort ensures that `v2` is treated as smaller than `v11`, whereas the default in `sort` (lexicographic sort) unintuitively treats `v2` as larger than `v11`.

	This is because a lexicographic sort does the sorting one character at a time, and does not treat the `2` and `11` as being a group of characters that represent a number. So since `1` comes before `2`, even version `v1000` is considered smaller than `v2` if you don't use natural sort.
</Alert>

<Exercise {criteria} {hints} />
