<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

Now, we have the tools to generate phylogenetic trees! Let's take a look at a slightly more challenging dataset: HIV-1 whole genome sequences collected from real people!

Let's go through the same sets of steps as before.

1. Run mutliple sequence alignment with MAFFT:
<Execute command="mafft hiv1_sequences.fas > hiv1_sequences.MSA.fas" inline />

2. Generate an unrooted tree with FastTree:
<Execute command="FastTree -gtr -nt -gamma hiv1_sequences.MSA.fas > tree_file.nwk
" inline />

3. Generate a rooted tree with LSD2:
<Execute command="lsd2 -i tree_file.nwk -d hiv1_dates.txt -l -1 -o lsd2_out" inline />

Now, let's take a look at the trees to ensure that they reflect what we would expect from the data.
