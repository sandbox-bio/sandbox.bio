<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

Sometimes we have a large collection of files and we actually just want to put them all together into a single file.

We can do this easily with the `bqtools cat` command.

<Execute command="bqtools cat fastq/sample*.vbq -o merged.vbq" />

Let's take a look at the number of reads in the merged file:

<Execute command="bqtools count merged.vbq" />

Great! They're all there.
