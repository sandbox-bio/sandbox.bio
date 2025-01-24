<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

**Generate a phylogenetic tree using FastTree**

Now that we have run multiple sequence alignment on our dataset, we can use [FastTree](https://morgannprice.github.io/fasttree/) to generate our phylogenetic tree.

1. Try <Execute command="FastTree" inline /> to
take a look at the usage instruction of FastTree.

- We will want to use a Generalized Time-Reversible (GTR) model which describes the rate of evolutionary charges between bases. So we will incorporate the flag '-gtr' in our command. 

- We also want to specify that our alignment data is made up of nucleotides, and not amino acids, so we will also include '-nt' in our command. 

2. Try <Execute command="FastTree -gtr -nt -gamma hiv1_aligned.txt > tree_file.nwk
" inline /> to generate our phylogenetic tree.
