<script>
// Solution:
//    seqkit grep --by-seq --pattern AACCGGUU hairpins.fa > pattern.fa

import Link from "$components/Link.svelte";
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>pattern.fa</code> exists",
	checks: [{
		type: "file",
		path: "pattern.fa",
		action: "exists"
	}]
},
{
	name: "File <code>pattern.fa</code> contains sequences from <code>hairpins.fa</code> that contain AACCGGUU",
	checks: [{
		type: "file",
		path: "pattern.fa",
		action: "contents",
        commandObserved: "echo 20",
        commandExpected: "cat pattern.fa | wc -l"  // much faster check than running seqkit seq
	}]
}];
</script>

Using <Link href="https://bioinf.shenwei.me/seqkit/usage/#grep">seqkit grep</Link>, create the file `pattern.fa`, which contains sequences from `hairpins.fa` that have the sequence `AACCGGUU`.

<Exercise {criteria} />
