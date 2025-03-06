<script>
import Link from "$components/Link.svelte";
</script>

Comparing genomes between different species or strains can reveal fascinating insights about evolution, structural variation, and genome organization. In this tutorial, we'll walk through the process of aligning two genomes and then visualizing their relationship using a circos plot in Circa.

We'll use data from <Link href="https://pmc.ncbi.nlm.nih.gov/articles/PMC10936905">Ament-Velásquez, et al.</Link>, a 2024 paper studying genome evolution in fungi. The paper compares genomes of <i>Podospora anserina</i> and related species, providing an excellent example of how genome alignments can illuminate evolutionary relationships. We will reproduce one of the circos plots in <Link href="https://pmc.ncbi.nlm.nih.gov/articles/PMC10936905/#evae034-F1">figure 1 (upper right corner)</Link> because it shows that <i>P. pseudocomata</i> aka CBS415.72 has an interesting set of alignments to <i>P. anserina</i>, with chromosome 1 and 2 from <i>P. anserina</i> each split into two pieces in <i>P. pseudocomata</i>.

We will align the two genomes to each other using MUMmer's `nucmer` tool, wrangle the output into a CSV file, and then visualize the results using Circa. Here is the visualization we will reproduce:

<Link href="https://circa.omgenomics.com/app/plot/gallery/aligned_genomes" />
