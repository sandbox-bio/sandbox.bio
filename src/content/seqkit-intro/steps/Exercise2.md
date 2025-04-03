<script>
// Solution:
//    seqkit seq --only-id --name NA12878.fastq > read_names.txt

import Link from "$components/Link.svelte";
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>read_names.txt</code> exists",
	checks: [{
		type: "file",
		path: "read_names.txt",
		action: "exists"
	}]
},
{
	name: "File <code>read_names.txt</code> contains SRR IDs of each read in <code>NA12878.fastq</code>",
	checks: [{
		type: "file",
		path: "read_names.txt",
		action: "contents",
        commandObserved: "echo 1000",
        commandExpected: "grep '^SRR' read_names.txt | grep -v 'E7VJJ' | wc -l"  // much faster check than running seqkit seq
	}]
}];
</script>

Using <Link href="https://bioinf.shenwei.me/seqkit/usage/#seq">seqkit seq</Link>, create the file `read_names.txt`, which contains the sequencing read names of the file `NA12878.fastq`. It should only contain the `SRR` IDs and not the full description.

<Exercise {criteria} />
