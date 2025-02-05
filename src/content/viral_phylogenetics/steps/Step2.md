<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

**Generate an unrooted phylogenetic tree using FastTree**

Now that we have run multiple sequence alignment on our dataset, we can use [FastTree](https://morgannprice.github.io/fasttree/) to generate our phylogenetic tree.

1. Try <Execute command="FastTree" inline /> to
take a look at the usage instruction of FastTree.

2. Now, try <Execute command="FastTree -gtr -nt -gamma hiv1_sequences.MSA.fas > tree_file.nwk
" inline /> to generate our phylogenetic tree.

Let's make some sense of this command:

- FastTree uses the [Jukes-Cantor](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)) or Generalized Time-Reversible ([GTR(https://en.wikipedia.org/wiki/Substitution_model#Generalised_time_reversible) models of nucleotide evolution. In our command, we will specifify that we want to use the generalized time-reversible model and incorporate the flag `-gtr` in our command. 

- We also want to specify that our alignment data is made up of nucleotides, and not amino acids, so we will also include `-nt` in our command.

- `-gamma` is an optional flag that allows for rescaling of the branch lengths and computation of a Gamma2-based liklihood

- `hiv1_sequences.MSA.fas > tree_file.nwk` tells FastTree to take in our multiple sequence alignment file (from Step 1) and output our phylogentic tree to a file called `tree_file.nwk`. A `nwk` file is a Newick format file, which is often used to represent phylogenetic trees. It is a text-based way to represent the tree structure.

