<script>
import Execute from "$components/Execute.svelte";
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>ids.fasta</code> exists",
	checks: [{
		type: "file",
		path: "ids.fasta",
		action: "exists"
	}]
},
{
	name: "File <code>ids.fasta</code> contains the FASTA sequences of sequences with IDs stored in <code>ids.txt</code>.",
	checks: [{
		type: "file",
		path: "ids.fasta",
		action: "contents",
		commandObserved: "wc ids.fasta",
        commandExpected: `echo '  21  198 2617 ids.fasta'`
	}]
}];
</script>

The file `ids.txt` contains a column of IDs identified as interesting in some way.

<Execute command="cat ids.txt" />

Use `blastdbcmd` to extract just those sequence records from the `orf_trans` database as a FASTA file named `ids.fasta`. Again, browsing the BLAST+ manual and the output of `blastdbcmd -help` will be useful.

<Exercise {criteria} />
