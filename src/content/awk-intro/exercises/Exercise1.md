<script>
// Solution:
//    awk -F "\t" '{ if($2 >= 5) print }' orders.tsv > large_orders.tsv

import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>large_orders.tsv</code> exists",
	checks: [{
		type: "file",
		path: "large_orders.tsv",
		action: "exists"
	}]
},
{
	name: "File <code>large_orders.tsv</code> contains orders where the same item was ordered >= 5 times (print all columns; no need to skip the header line)",
	checks: [{
		type: "file",
		path: "large_orders.tsv",
		action: "contents",
		commandExpected: `awk -F "\\t" '{ if($2 >= 5) print }' orders.tsv`
	}]
}];
</script>

Create a tab-separated file called `large_orders.tsv` that contains all the lines from the file `orders.tsv` where someone ordered the same item 5 times or more (note that column 2 contains the item count).

<Exercise {criteria} />
