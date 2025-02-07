<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

**Generate an unrooted phylogenetic tree using FastTree**

Now that we have run multiple sequence alignment on our dataset, we can use [FastTree](https://morgannprice.github.io/fasttree/) to generate our phylogenetic tree.

1. Try <Execute command="FastTree" inline /> to
take a look at the usage instructions.

2. Now, try <Execute command="FastTree -gtr -nt -gamma hiv1_sequences.MSA.fas > tree_file.nwk
" inline /> to generate our phylogenetic tree.

Let's make some sense of this command:

- `-gtr` indicates we will be using a [Generalized Time-Reversible](https://en.wikipedia.org/wiki/Substitution_model#Generalised_time_reversible) (GTR) model of evolution for our tree. FastTree can be run with either the [Jukes-Cantor](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)) or GTR models.

- `-nt` specifies that our alignment data is made up of nucleotides, and not amino acids

- `-gamma` is an optional flag that allows for rescaling of the branch lengths and computation of a Gamma2-based liklihood

- `hiv1_sequences.MSA.fas > tree_file.nwk` tells FastTree to take in our multiple sequence alignment file (from Step 1) and output our phylogentic tree to a file called `tree_file.nwk`. A `.nwk` file is a Newick format file, which is often used to represent phylogenetic trees. It is a text-based way to represent the tree structure

After running the above command, we will have our unrooted phylogenetic tree. Use <Execute command="head -10 tree_file.nwk" inline /> to view the first 10 lines of the tree file. You can read more about Newick file formats [here](https://en.wikipedia.org/wiki/Newick_format).

Now, use <Execute command="downlaod tree_file.nwk" inline /> to store the file locally. Navigate to [Taxonium](https://taxonium.org/?xType=x_dist) to upload the file, and view your unrooted phylogenetic tree. 
