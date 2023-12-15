<script>
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
</script>

<Alert>
	This tutorial is an interactive version of the <Link href="http://quinlanlab.org/tutorials/bedtools.html">bedtools tutorial</Link> developed by the <Link href="http://quinlanlab.org/">Quinlan Lab</Link>. The contents are the same, but the data was subsampled so it can be analyzed in your browser.
</Alert>

Our goal is to work through examples that demonstrate how to explore, process and manipulate genomic interval files (e.g., `BED`, `VCF`, `BAM`) with the `bedtools` software package.

Some of our analysis will be based upon the <Link href="https://science.sciencemag.org/content/337/6099/1190">Maurano et al.</Link> exploration of DnaseI hypersensitivity sites in hundreds of primary tissue types.

This tutorial is merely meant as an introduction to whet your appetite. There are many, many more tools and options than presented here. We therefore encourage you to read the <Link href="https://bedtools.readthedocs.io/en/latest/">bedtools documentation</Link>.
