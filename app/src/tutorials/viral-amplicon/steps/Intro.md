<script>
import { onMount } from "svelte";
import { CLI } from "./terminal/cli";
import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";
import Execute from "./components/Execute.svelte";
</script>

<Alert>
	**Note**: Here we assume you completed the <Link href="/tutorials?id=bowtie2-intro">bowtie2</Link> and <Link href="/tutorials?id=samtools-intro">samtools</Link> tutorials. The analysis workflow we will be using is adapted from <Link href="https://github.com/ucsd-ccbb/C-VIEW">UC San Diego's C-VIEW pipeline</Link>
</Alert>

During the COVID-19 panemic, obtaining viral genome sequences from samples collected from patients or the environment (e.g. wastewater) enabled public health officials to track the evolution of the SARS-CoV-2 over time, identify mutations of interest/concern, and detect the prevalence and growth of new viral lineages as the pandemic progressed. Because a high-confidence reference genome sequence was obtained very early on in the pandemic, researchers were able to perform "amplicon sequencing": a highly-targeted approach for sequencing specific pre-known regions of interest across the viral genome.

In this tutorial, we'll explore how to call variants and reconstruct a "consensus" genome sequence from a given SARS-CoV-2 amplicon sequencing dataset.