<script>
// Solution:
//    bedtools flank -l 2 -r 2 -i exons.bed -g genome.txt > splice-sites.bed

import Exercise from "$components/Exercise.svelte";
import Alert from "$components/Alert.svelte";
import Link from "$components/Link.svelte";

let criteria = [
{
	name: "File <code>splice-sites.bed</code> contains a list of splice sites",
	checks: [{
		type: "file",
		path: "splice-sites.bed",
		action: "contents",
		commandExpected: "bedtools flank -l 2 -r 2 -i exons.bed -g genome.txt"
	}]
}
];
</script>

Create intervals representing the canonical 2bp splice sites on either side of each exon (don't worry about excluding splice sites at the first or last exon). Use the files `exons.bed` and `genome.txt` as input, and output your result to `splice-sites.bed`.

<Alert>
	**Hint**: have a look at the <Link href="https://bedtools.readthedocs.io/en/latest/content/tools/flank.html">`bedtools flank`</Link> tool.
</Alert>

<Exercise {criteria} />
