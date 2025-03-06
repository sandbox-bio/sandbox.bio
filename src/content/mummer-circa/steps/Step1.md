<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

We preloaded data from <Link href="https://pmc.ncbi.nlm.nih.gov/articles/PMC10936905">Ament-Velásquez, et al.</Link>. We got this data from the <Link href="https://datadryad.org/stash/dataset/doi:10.5061/dryad.1vhhmgr0j">supplementary data</Link> by downloading the big `Dryad.zip` file and finding the two files we needed inside the "Assemblies" folder.

To show the files we are working with, click on or type <Execute command={"ls"} inline /> into the command line.

It's always good to take a look at both files to confirm they look like proper FASTA files before we start.

<Execute command={"head Podan2_AssemblyScaffoldsmt.fa"} />

<Execute command={"head CBS415.72m.nice_mt.fa"} />

Let's start by aligning two genomes using the MUMmer4 package. We use `nucmer` because that is MUMmer's aligner for nucleotide sequences ("nuc"), i.e. DNA or RNA sequences, and here we are aligning DNA.

With nucmer, we need to pick a "reference" and a "query". Usually the reference is the more complete or established genome, while the query is the new one. In this paper, they are comparing many different genomes all to Podan2 (an existing genome assembly that predates the paper), so in this case Podan2 is the clear pick for the reference, which makes the query CBS415.72.

You can see Nucmer's expected syntax by simply trying to run it (i.e. <Execute command={"nucmer"} inline />), but just note that the reference should come before the query, i.e. `nucmer <options> <reference.fa> <query.fa>`.

<Execute command={"time nucmer -maxmatch -l 100 -c 1000 \\ --prefix CBS415_to_Podan \\ Podan2_AssemblyScaffoldsmt.fa \\ CBS415.72m.nice_mt.fa"} />

That should take a few minutes, so grab a cup of coffee, then come back and view the output:

<Execute command={"head CBS415_to_Podan.delta"} />

That's not very human-readable, but we can use MUMmer's utilities to view it:
<Execute command={"show-coords CBS415_to_Podan.delta | head"} />
