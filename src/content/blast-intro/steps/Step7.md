<script>
import Link from "$components/Link.svelte";
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>p450s_blastp.txt</code> exists",
	checks: [{
		type: "file",
		path: "p450s_blastp.txt",
		action: "exists"
	}]
},
{
	name: "File <code>p450s_blastp.txt</code> is the equivalent of <code>p450s_blastp.blast</code>, but using output format 6",
	checks: [{
		type: "file",
		path: "p450s_blastp.txt",
		action: "contents",
		commandObserved: "grep 'CP1A1_CANFA' p450s_blastp.txt",
        commandExpected: `echo -e "sp|P56590|CP1A1_CANFA\\tYHR007C\\t24.101\\t278\\t174\\t9\\t240\\t492\\t236\\t501\\t9.36e-17\\t74.7"`
	}]
}];
</script>

Use the `blast_formatter` tool to convert the output format 11 file above into an output format 6 called `p450s_blastp.txt`. You may find browsing the <Link href="https://www.ncbi.nlm.nih.gov/books/NBK569843/">NCBI BLAST+ manual</Link> and the output of `blast_formatter -help` to be informative.

<Exercise {criteria} />
