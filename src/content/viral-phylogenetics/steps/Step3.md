<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

Now that we have run MSA on our dataset, we can perform Phylogenetic Inference: using the MSA and mathematical models of evolution, we can identify the evolutionary tree (i.e. phylogeny) that best describes the evolutionary relationships between the sequences in our MSA. It is important to note that the phylogeny we produce in this step will be unrooted. This means that while we can determine the relative evolutionary relationships between sequences, we cannot "root" the tree (and thus identify the date of the MRCA) because the directionality of time is not known. In an unrooted tree, branch lengths are in units of mutations instead of in units of time.

We will use [FastTree](https://morgannprice.github.io/fasttree/) to generate an unrooted phylogenetic tree to assess relative ancestral relationships between our samples.

1. Try <Execute command="FastTree" inline /> to take a look at the usage instructions.

2. Now, generate our phylogenetic tree: 

<Execute command="FastTree -nt ViralMSA_Out/sarscov2_sequences.fas.aln > sarscov2_sequences.unrooted_tree.nwk" />

Let's make some sense of this command:

- `-nt` specifies that our alignment is of **n**ucleo**t**ides and not amino acids

- `sarscov2_sequences.msa.fas > sarscov2_sequences.unrooted_tree.nwk` tells FastTree to take in our multiple sequence alignment file (from Step 1) as input and to output the unrooted phylogenetic tree to a file called `sarscov2_sequences.unrooted_tree.nwk` in the same directory. A `.nwk` file is in Newick format, which is often used to represent phylogenetic trees. It is a text-based way to represent the tree structure. You can read more about Newick format [here](https://en.wikipedia.org/wiki/Newick_format).

**Optional Options:**

Note: The following options may be included in the command to improve accuracy at the cost of increased runtime.

- `-gtr` implements use of a [Generalized Time-Reversible](https://en.wikipedia.org/wiki/Substitution_model#Generalised_time_reversible) (GTR) model of evolution for our tree. FastTree can be run with either the [Jukes-Cantor](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)) or GTR model. 

- `-gamma` allows for rescaling of the branch lengths and computation of a Gamma2-based likelihood


3. Use `head` to view the first 10 lines of the tree file:

<Execute command="head -10 sarscov2_sequences.unrooted_tree.nwk" /> 

4. Now, let's quickly visualize how this information makes a tree in the terminal using `nw_display`:

<Execute command="nw_display sarscov2_sequences.unrooted_tree.nwk" />

Why might we want to create an unrooted tree, over a rooted tree (with a common ancestor)?

<Quiz
	id="step2-quiz1"
	choices={[
		{ valid: false, value: `To determine the evolutionary direction and ancestral lineage of species` },
		{ valid: true, value: `To analyze relationships without assuming a common ancestor or direction of evolution` },
		{ valid: false, value: `Unrooted trees are always more accurate than rooted trees` },
		{ valid: false, value: `To define the exact point in time when species diverged` },
    ]}>
	<span slot="prompt"></span>
</Quiz>

5. We can alternatively visualize our tree using webtools. First, use `download` to store the file locally:

<Execute command="download sarscov2_sequences.unrooted_tree.nwk" /> 

6. Navigate to [Taxonium](https://taxonium.org/?xType=x_dist) to upload the file, and view your unrooted phylogenetic tree. 
