<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

In the previous step we created an unrooted phylogenetic tree through Phylogenetic Inference. However, because we have access to the collection of dates of the sequences in our SARS-CoV-2 dataset, we can "root" the tree (find the most likely position of the MRCA) and "date" the tree (scale the branch lengths to be in units of time). We will use [LSD2](https://github.com/tothuhien/lsd2) to generate a rooted tree, and we will use an outgroup to help us do this. Known organisms that are distantly related to the species of interest can act as outgroups (i.e. references) when inferring a rooted tree, which can help us perform more accurate rooting and dating. In our case, we will use a RaTG13 bat coronavirus sequence as our outgroup. 

1. Try <Execute command="LSD2 -help" inline /> to take a look at the usage instruction of LSD2.

2. Now, to generate our rooted tree, use <Execute command="lsd2 -i sarscov2_sequences.unrooted_tree.nwk -d sarscov2_dates.txt -g sarscov2_outgroup.txt -G -l -1 -o lsd2_out" inline /> 

The above command incorporates the following flags:

- `-i` specifies the input file, which is our unrooted phylogenetic tree from Step 2
- `-d` specifies the file with sequences dates, which is essential for rooting
- `-g` specifies the file with outgroup sequences
- `-G` removes the outgroups from the tree (uses it to root, but does not show it on the tree)
- `-o` specifies the name of our output file

3. Now, we have a rooted tree stored in a file called `lsd2_out.nwk`. Like in Step 2, we can view the first 10 lines of the Newick file at the command line with <Execute command="head -10 lsd2_out.nwk" inline />

4. Let's visualize in the terminal using <Execute command="nw_display - < lsd2_out.nwk" inline />

5. Again, we can download the file with <Execute command="download lsd2_out.nwk" inline /> so that it can be uploaded and better visualized in [Taxonium](https://taxonium.org/?xType=x_dist). 

Take a look at the `lsd2_out` log file. When did the MRCA exist?

<Quiz
	id="step3-quiz1"
	choices={[
		{ valid: false, value: `2019-08-12` },
		{ valid: false, value: `2019-03-11` },
		{ valid: true, value: `2020-03-02` },
		{ valid: false, value: `2020-01-28` },
    ]}>
	<span slot="prompt"></span>
</Quiz>

WHO declared COVID-19 a pandemic on March 11, 2020. Does our MRCA date to before or after this day?

<Quiz
	id="step3-quiz2"
	choices={[
		{ valid: true, value: `before` },
		{ valid: false, value: `after` },
    ]}>
	<span slot="prompt"></span>
</Quiz>

We used 10 SARS-CoV-2 sequences to generate this rooted tree. Which statement is true about this approach?

<Quiz
	id="step3-quiz3"
	choices={[
		{ valid: false, value: `It guarantees an accurate estimate of the virus's mutation rate.` },
		{ valid: false, value: `It captures the full geographic and temporal diversity of the virus.` },
		{ valid: true, value: `It may result in a tree that does not accurately reflect the true evolutionary history.` },
		{ valid: false, value: `It eliminates the need for an outgroup to root the tree.` },
    ]}>
	<span slot="prompt"></span>
</Quiz>

The fact that the date of the MRCA is in April 2019 reflects the fact that SARS-CoV-2 is [incredibly difficult to root](https://doi.org/10.1126/science.abp8337). Outgroup rooting does a decent job getting us close to the correct MRCA, but more sophisticated methods are needed to improve the accuracy of SARS-CoV-2 MRCA. 