<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

First, let's look at our input files:

<Execute command="ls" />

We have two FASTA files:

- H_pylori26695_Eslice.fasta: The reference strain
- H_pyloriJ99_Eslice.fasta: The query strain

These FASTA files come from the <Link href="https://mummer4.github.io/tutorial/tutorial.html">MUMmer tutorial</Link>.

Let's set up some variables to make our commands cleaner:

<Execute command={`REF=H_pylori26695_Eslice.fasta
QUERY=H_pyloriJ99_Eslice.fasta
NAME=\${REF%.fasta}_to_\${QUERY%.fasta}`} />

Let's start by aligning two genomes using the MUMmer4 package. We use `nucmer` because that is MUMmer's aligner for nucleotide sequences ("nuc"), i.e. DNA or RNA sequences, and here we are aligning DNA.

<Execute command="nucmer --prefix $NAME $REF $QUERY" />
(This might take a minute)

This will create a delta file containing the alignments. The `--prefix` option sets the prefix for output files.

Let's look at the first few lines of the delta file:

<Execute command="head $NAME.delta" />

The delta file format is not very human-readable, so in the next step we'll convert it to a more useful format.
