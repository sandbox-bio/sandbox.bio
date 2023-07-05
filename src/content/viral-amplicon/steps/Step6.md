<script>
import Execute from "$components/Execute.svelte";
</script>

Now that we have a pile-up file, we'll use `ivar` to call variants: all positions in which our sample deviates from the reference genome. Run the following:

<Execute command="cat pileup.txt | \ ivar variants \ -r $REF_FASTA \ -g $REF_GFF \ -p variants.tsv -m 10" />

Let's break this seemingly complex command into its individual components to make some sense of it:

- First, we're calling the following command: `cat pileup.txt`
  - `cat` prints the contents of `pileup.txt` to standard output
- Then, we're using the pipe character `|` to pipe the output of `cat` (via standard output) to `ivar` (via standard input)
- The `ivar` command we're piping the pile-up stream into is the following: `ivar variants -r $REF_FASTA -g $REF_GFF -p variants.tsv -m 10`
  - `variants` tells `ivar` that we want to use its variant-calling functionality
  - `-r $REF_FASTA` tells `ivar` that we want to use the file `$REF_FASTA` as our reference genome FASTA
  - `-g $REF_GFF` tells `ivar` that we want to use the file `$REF_GFF` as our reference genome's annotation GFF
  - `-p variants.tsv` tells `ivar` that we want to write the output variants to the file `variants.tsv`
    - The output file is a TSV file: **T**ab-**S**eparated **V**alues
  - `-m 10` tells `ivar` that we want to only output variants that have a **m**inimum depth of 10
    - In other words, variants in which there are at least 10 reads that cover that position of the reference genome

After running the above command, we will have sucessfully called variants and written the results to the file `variants.tsv`.

To see the first few lines of the variants output file, run the following:

<Execute command="head -n 5 variants.tsv" />

This variants TSV file is one of the key results of our amplicon sequencing analysis.
