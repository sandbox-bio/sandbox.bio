<script>
// Solution:
//    awk -F "\t" '{ if($3 == "Bottled Water" && $2 > largest) largest=$2 } END { print largest }' orders.tsv > largest_order.tsv

import Exercise from "components/Exercise.svelte";

let criteria = [
{
	name: "File <code>largest_order.tsv</code> exists",
	checks: [{
		type: "file",
		path: "largest_order.tsv",
		action: "exists"
	}]
},
{
	name: "File <code>largest_order.tsv</code> contains the largest order of <code>Bottled Water</code> by a single customer",
	checks: [{
		type: "file",
		path: "largest_order.tsv",
		action: "contents",
		commandExpected: `awk -F "\t" '{ if($3 == "Bottled Water" && $2 > largest) largest=$2 } END { print largest }' orders.tsv`
	}]
}];
</script>

Create a file called `large_order.tsv` that outputs the largest order of `Bottled Water` by a single customer. You only need to output the number of water bottles, so the file should only contain an integer. Note that column 3 contains the item name, and column 2 contains the quantity of each item.

<Exercise {criteria} />
