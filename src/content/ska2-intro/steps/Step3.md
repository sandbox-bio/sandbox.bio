<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

If no further pruning is needed, we can just use the alignment from the previous step to build a phylogeny. There are various tools for doing this, such as `VeryFastTree` and `ggtree`. We will use a very simple Python script that makes use of the <Link href="https://biopython.org/">BioPython library</Link>, that runs quickly in the browser. It will read the alignment, create a tree, and print it in terminal. Run it executing the following:

<Execute command="python3 create_tree.py output/ska_alignment.aln" />

This shows that some of the assemblies are more closely related than the others.
