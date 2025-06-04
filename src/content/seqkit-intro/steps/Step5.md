<script>
import Execute from "$components/Execute.svelte";
</script>

If your file contains duplicate sequences, you can use `rmdup` to remove them:

<Execute command="seqkit rmdup hairpins.fa > deduplicated.fa" />

To see the duplicate sequences:

<Execute command="diff deduplicated.fa hairpins.fa" />

> 💡 You can use the `--by-name` flag to check for duplicates using the full sequence name instead of just the ID (e.g. `cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop` instead of just the ID `cel-let-7`).
>
> You can also use `--by-seq` to check for duplicates using only the sequences, regardless of whether the IDs are different.
