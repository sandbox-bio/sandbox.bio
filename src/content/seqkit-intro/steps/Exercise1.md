<script>
// Solution:
//    seqkit seq -w 80 hairpins.fa > formatted.fa

import Link from "$components/Link.svelte";
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>formatted.fa</code> exists",
	checks: [{
		type: "file",
		path: "formatted.fa",
		action: "exists"
	}]
},
{
	name: "File <code>formatted.fa</code> has 80 bases per line",
	checks: [{
		type: "file",
		path: "formatted.fa",
		action: "contents",
        commandObserved: "echo 81",
        commandExpected: "sed '2q;d' formatted.fa | wc -c"  // much faster check than running seqkit seq
	}]
}];
</script>

Using <Link href="https://bioinf.shenwei.me/seqkit/usage/#seq">seqkit seq</Link>, create the file `formatted.fa`, which formats `hairpins.fa` to show 80 bases per line instead of the current 60 bases per line.

<Exercise {criteria} />
