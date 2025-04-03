<script>
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>p450s_blastp.blast</code> exists",
	checks: [{
		type: "file",
		path: "p450s_blastp.blast",
		action: "exists"
	}]
},
{
	name: "File <code>p450s_blastp.blast</code> contains the output of aligning sequences from <code>p450s.fasta</code> to the <code>orf_trans</code> database.",
	checks: [{
		type: "file",
		path: "p450s_blastp.blast",
		action: "contents",
		commandObserved: "grep 'MSTSAMELLLTATIFCLVLWVVRIFRPQVPKGLKSPPGPWGWPLIG' p450s_blastp.blast",
        commandExpected: `echo '              seq-data ncbieaa "MSTSAMELLLTATIFCLVLWVVRIFRPQVPKGLKSPPGPWGWPLIG'`
	}]
}];
</script>

Align the file `p450s.fasta` to the database we created in previous steps.

When doing the search, use an E-value cutoff of `1e-6`, keep the top one target sequences, and produce an output file called `p450s_blastp.blast` in output format 11 (`-outfmt 11`).

<Exercise {criteria} />
