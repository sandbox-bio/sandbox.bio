<script>
import Link from "$components/Link.svelte";
</script>

In this tutorial, we'll learn how to align and compare two bacterial genomes using MUMmer. We'll use two different strains of _Helicobacter pylori_ - the reference strain 26695 and strain J99.

_H. pylori_ is a bacterium that colonizes the human stomach and can cause ulcers and gastric cancer. By comparing different strains, we can better understand how this pathogen evolves and adapts to different hosts.

We'll use MUMmer's `nucmer` tool to align the genomes, process the output into a readable format, and prepare the data for visualization. The workflow will be:

1. Align genomes with `nucmer`
2. Convert alignment to coordinates with `show-coords`
3. Process the data into CSV format
4. Generate reference files that will set up the coordinate system for visualization
5. Use Circa to make a circos plot

Let's get started!
