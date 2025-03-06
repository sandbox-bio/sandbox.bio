<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

In the previous step we created an unrooted phylogenetic tree, now, we will use [LSD2](https://github.com/tothuhien/lsd2) to create a rooted tree, and use an outgroup to help us do this. Known organisms that are distantly related to the species of interest, can act as outgroups or references when building a rooted tree. 

1. Try <Execute command="LSD2 -help" inline /> to take a look at the usage instruction of LSD2.

2. Now, to generate our rooted tree, use <Execute command="lsd2 -i sarscov2_sequences.unrooted_tree.nwk -d sarscov2_dates_lsd2.txt -g sarscov2_outgroup.txt -G -l -1 -o lsd2_out" inline /> 

The above command incorporates the following flags:

- `-i` specifies the input file, which is our unrooted phylogenetic tree from Step 2
- `-d` specifies the file with sequences dates, which is essential for rooting
- `-g` specifies the file with outgroup sequences
- `-G` removes the outgroups from the tree (uses it to root, but does not show it on the tree)
- `-o` specifies the name of our output file

3. Now, we have a rooted tree stored in a file called `lsd2_out.tree.result.nwk`. Like in Step 2, we can view the first 10 lines of the Newick file at the command line with <Execute command="head -10 phylogenetic.tree.result.nwk" inline />

4. Let's visualize in the terminal using <Execute command="nw_display - < phylogenetic.tree.result.nwk" inline />

5. Again, we can download the file with <Execute command="download phylogenetic.tree.result.nwk" inline /> so that it can be uploaded and better visualized in [Taxonium](https://taxonium.org/?xType=x_dist). 
