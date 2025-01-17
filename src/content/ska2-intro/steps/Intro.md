<script>
import Link from "$components/Link.svelte";
import Image from "$components/Image.svelte";
</script>

<Link href="https://github.com/bacpop/ska.rust">SKA</Link> is a tool for comparing small and highly similar genomes using <Link href="https://www.biorxiv.org/content/early/2018/10/25/453142">split _k_-mers</Link>.

This tutorial (adapted from [this guide](https://www.bacpop.org/guides/building_trees_with_ska)) will explain how to use SKA to build a phylogenetic tree for some _Escherichia coli_ lineages in a few minutes. Although SKA is tailored more towards analysing variation within a lineage, tree-building ends up working fine for the whole species but requires more memory.

<Image src="/data/ska2-intro/ska_logo.png" />
