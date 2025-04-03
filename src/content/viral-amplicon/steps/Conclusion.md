<script>
import Link from "$components/Link.svelte";
</script>

Hopefully this tutorial gave you some insights into how viral amplicon sequence data are analyzed! This is the exact general workflow that was used to produce over 16 million SARS-CoV-2 genome sequences across the world, which were immensely helpful for tracking the spread of COVID-19.

We happened to use the ViralConsensus pipeline produced by <Link href="https://niema.net/">Niema Moshiri at UC San Diego</Link>, but in practice, you may want to see how swapping out tools for each of these steps change the runtime and results of the analysis. If you're curious in learning about how different stages of this pipeline could be modified, check out <Link href="https://doi.org/10.1038/s41598-022-09035-w">Moshiri et al., Scientific Reports 2022</Link>.
