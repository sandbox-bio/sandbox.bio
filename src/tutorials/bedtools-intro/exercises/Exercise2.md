<script>
// Solution:
//    bedtools intersect -a <(bedtools makewindows -g genome.txt -w 500000) -b exons.bed -c > windows.exons.bedg

import Alert from "components/Alert.svelte";
import Exercise from "components/Exercise.svelte";

let criteria = [
{
	name: "File <code>windows.bed</code> contains a list of all regions of 500kb in the genome",
	checks: [{
		type: "file",
		path: "windows.bed",
		action: "contents",
		commandExpected: "bedtools makewindows -g genome.txt -w 500000",
	}]
},
{
	name: "File <code>windows.exons.bedg</code> contains a list of each 500kb interval and how many exons were found within that region",
	checks: [{
		type: "file",
		path: "windows.exons.bedg",
		action: "contents",
		commandExpected: "bedtools intersect -a <(bedtools makewindows -g genome.txt -w 500000) -b exons.bed -c"
	}]
}

];
</script>

Count how many exons occur in each 500kb interval ("window") in the human genome. Use the files `exons.bed` and `genome.txt` as input.

<Alert>
	**Hint**: have a look at the `bedtools makewindows` tool.
</Alert>

<Exercise {criteria} />
