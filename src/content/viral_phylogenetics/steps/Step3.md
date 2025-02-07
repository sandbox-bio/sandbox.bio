<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

In the previous step we created an unrooted phylogenetic tree, now, we will use [LSD2](https://github.com/tothuhien/lsd2) to create a rooted tree. LSD2 is a program that uses least-squares methods to estimate rates and dates.  

1. Try <Execute command="LSD2 -help" inline /> to take a look at the usage instruction of LSD2.

2. Now, to generate our rooted tree, use <Execute command="lsd2 -i tree_file.nwk -d hiv1_dates.txt -r a -l -1 -u 0 -q 0.2 -R 365 -t 0.00000000010000000000 -v 1 -s 9182" inline /> 

The above command incorporates the following flags:

- `-i` specifies the input file, which is our unrooted phylogenetic tree from Step 2
- `-d` specifies the file with sequences dates, which is essential for rooting
- `-r a` sets the rooting method to automatic
- `-l -1` sets the log-liklihood optimization rate, -1 is the default
- `-u 0` defined the upper bound for optimization rate
- `-q 0.2` sets the quantile threshold for outlier removal
- `-R 365` specifies that the dates given are in years (365 days/year)
- `-t 0.00000000010000000000` defines the initial time estimate or starting molecular clock rate
- `-v 1` sets the verbosity level to 1
- `-s 9182` sets the random seed for reproducibility

3. Now, we have a rooted tree stored in a file called `phylogenetic.tree.result.nwk`. Like in Step 2, we can view the first 10 lines of the Newick file at the command line with <Execute command="head -10 phylogenetic.tree.result.nwk" inline />

4. Let's visualize in the terminal using <Execute command="nw_display - < phylogenetic.tree.result.nwk" inline />

5. Again, we can download the file with <Execute command="download phylogenetic.tree.result.nwk" inline /> so that it can be uploaded and better visualized in [Taxonium](https://taxonium.org/?xType=x_dist). 
