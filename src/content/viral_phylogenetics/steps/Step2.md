<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

Now that we have run multiple sequence alignment on our dataset, we can use [FastTree](https://morgannprice.github.io/fasttree/) to generate an unrooted phylogenetic tree to assess relative ancestral relationships between our samples.

1. Try <Execute command="FastTree" inline /> to
take a look at the usage instructions.

2. Now, try <Execute command="FastTree -gtr -nt -gamma sarscov2_sequences.msa.fas > sarscov2_sequences.unrooted_tree.nwk" inline /> to generate our phylogenetic tree.

Let's make some sense of this command:

- `-gtr` indicates we will be using a [Generalized Time-Reversible](https://en.wikipedia.org/wiki/Substitution_model#Generalised_time_reversible) (GTR) model of evolution for our tree. FastTree can be run with either the [Jukes-Cantor](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)) or GTR model.

- `-nt` specifies that our alignment data is made up of nucleotides, and not amino acids

- `-gamma` is an optional flag that allows for rescaling of the branch lengths and computation of a Gamma2-based liklihood

- `cov2_sequences.msa.fas > sarscov2_sequences.unrooted_tree.nwk` tells FastTree to take in our multiple sequence alignment file (from Step 1) and output our phylogentic tree to a file called `tree_file.nwk`. A `.nwk` file is a Newick format file, which is often used to represent phylogenetic trees. It is a text-based way to represent the tree structure.

3. After running the above command, we will have our unrooted phylogenetic tree. Use <Execute command="head -10 tree_file.nwk" inline /> to view the first 10 lines of the tree file. You can read more about Newick file formats [here](https://en.wikipedia.org/wiki/Newick_format).

4. Now, let's quickly visualize how this information makes a tree in the terminal using <Execute command="nw_display - < tree_file.nwk" inline />


**Why might be want to create an unrooted tree, over a rooted tree (with a common ancestor)?**

<Quiz
	id="step2-quiz1"
	choices={[
		{ valid: false, value: `To determine the evolutionary direction and ancestral lineage of species` },
		{ valid: true, value: `To analyze relationships without assuming a common ancestor or direction of evolution` },
		{ valid: false, value: `Because unrooted trees are always more accurate than rooted trees` },
		{ valid: false, value: `To define the exact point in time when species diverged` },
    ]}>
	<span slot="prompt"></span>
</Quiz>

5. We can visualize our tool even better using webtools. First, use <Execute command="download tree_file.nwk" inline /> to store the file locally.
   
6. Navigate to [Taxonium](https://taxonium.org/?xType=x_dist) to upload the file, and view your unrooted phylogenetic tree. 
