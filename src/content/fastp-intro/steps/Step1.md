<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

In this tutorial, we'll use a small subset of a <Link href="https://42basepairs.com/browse/gs/deepvariant/benchmarking/fastq/wgs_pcr_free/30x?file=HG004.novaseq.pcr-free.30x.R1.fastq.gz">publicly available dataset</Link> provided by the <Link href="https://github.com/google/deepvariant">DeepVariant</Link> team.

This data was obtained from a paired-end sequencing experiment of the <Link href="https://catalog.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=NA24143&Product=DNA">NA24143</Link> sample.

Use <Execute inline command="ls" /> to list the 2 data files.

In the next step, we'll run `fastp` to generate a basic QC report.
