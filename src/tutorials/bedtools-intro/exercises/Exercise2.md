<script>
// Solution:
//    bedtools makewindows -g genome.txt -w 500000 > windows.bed
//    bedtools intersect -a windows.bed -b exons.bed -c > windows.exons.bedg

import Alert from "../../Alert.svelte";
import Exercise from "../../Exercise.svelte";

let criteria = [
{
	name: "File <code>windows.bed</code> exists",
	checks: [{
		type: "file",
		path: "windows.bed",
		action: "exists"
	}]
},
{
	name: "File <code>windows.bed</code> contains a list of all regions of 500kb in the genome",
	checks: [{
		type: "file",
		path: "windows.bed",
		action: "contents",
		equal: "bedtools makewindows -g genome.txt -w 500000",
		output: "/shared/tmp/exercise2-windows.bed"
	}]
},
{
	name: "File <code>windows.exons.bedg</code> exists",
	checks: [{
		type: "file",
		path: "windows.exons.bedg",
		action: "exists"
	}]
},
{
	name: "File <code>windows.exons.bedg</code> contains a list of each 500kb interval and how many exons were found within that region",
	checks: [{
		type: "file",
		path: "windows.exons.bedg",
		action: "contents",
		equal: "bedtools intersect -a /shared/tmp/exercise2-windows.bed -b exons.bed -c"
	}]
}

];
</script>

Count how many exons occur in each 500kb interval ("window") in the human genome. Use the files `exons.bed` and `genome.txt` as input.

<Alert>
	**Hint**: have a look at the `bedtools makewindows` tool.
</Alert>

<Exercise {criteria} />
