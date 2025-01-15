<script>
import Link from "$components/Link.svelte";
</script>

The basic approach to building a tree with SKA is to generate a SNP alignment using split _k_-mers and then feed that to a tree building algorithm of choice.

Since SKA does not require a specifying a reference sequence to determine the SNPs, SKA gets around potential biases introduced by reference choice and is thus particularly well suited to analysing microbial genomes derived from outbreaks or pathogen surveillance.

In this tutorial SKA v0.3.11 is used. The data we are going to use are four examples of _E. coli_ hybrid Nanopore+Illumina assemblies isolated from travellers visiting the city of Vientiane in Laos. This is a subset from a larger dataset, available <Link href="https://zenodo.org/record/8172518/files/building_trees_with_ska.tar">in this Zenodo entry</Link>: SKA is able to work with much more data in a normal computer, but we will use those four data entries as example for this tutorial.
