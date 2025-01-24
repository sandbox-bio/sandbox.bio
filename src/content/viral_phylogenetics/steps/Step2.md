<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

**Generate our phylogenetic tree using FastTree**

Now that we have run multiple sequence alignment on our dataset, we can use [FastTree](https://morgannprice.github.io/fasttree/) to generate our phylogenetic tree.

Try <Execute command="FastTree --h" inline /> to
take a look at all the usage instruction of FastTree.

