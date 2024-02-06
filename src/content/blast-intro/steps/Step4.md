<script>
import Link from "$components/Link.svelte";
import Image from "$components/Image.svelte";
import Execute from "$components/Execute.svelte";
</script>

Now we need to determine what options we will use for the `blastp`. In particular, do we want to limit the number of HSPs and target sequences reported for each query? Because we're mostly interested in determining which proteins match others, we probably only need to keep one hit. **But each protein's best hit will likely be to itself**! So we'd better keep the top two with `-max_target_seqs 2` and only the best HSP per hit with `-max_hsps 1`.

For the output, we'll create a tab-separated output with comment lines (`-outfmt 7`) called `yeast_blastp_yeast_top2.txt`:

<Execute command={`blastp \\ -query orf_trans.fasta \\ -db orf_trans \\ -max_target_seqs 2 \\ -max_hsps 1 \\ -evalue 1e-6 \\ -outfmt '7 qseqid sseqid length qlen slen qstart qend sstart send evalue' \\ -out yeast_blastp_yeast_top2.txt`} />

The coded names&mdash;`qseqid`, `sseqid`, `length`, etc.&mdash;can be found by running `blastp -help`.
