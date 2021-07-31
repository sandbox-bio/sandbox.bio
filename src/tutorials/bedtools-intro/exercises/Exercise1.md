<script>
// Solution:
//    bedtools complement -i exons.bed -g genome.txt > notexons.bed

import Exercise from "../../Exercise.svelte";

let criteria = [
{
	name: "File <code>notexons.bed</code> exists",
	checks: [{
		type: "file",
		path: "notexons.bed",
		action: "exists"
	}]
},
{
	name: "File <code>notexons.bed</code> contains non-exonic regions",
	checks: [{
		type: "file",
		path: "notexons.bed",
		action: "contents",
		equal: "bedtools complement -i exons.bed -g genome.txt"
	}]
}];
</script>

Create a BED file called `notexons.bed` that contains all of the intervals in the genome that are NOT exonic. Use the files `exons.bed` and `genome.txt` as input.

<Exercise {criteria} />
