<script>
import Link from "components/Link.svelte";
import Execute from "components/Execute.svelte";
</script>

**1. The sequencing reads**

The reads have been preloaded to your working directory. We will be analyzing paired-end sequence data, so we have two FASTQ files: `reads_R1.fq.gz` (representing the "Read 1" reads of each read-pair) and `reads_R2.fq.gz` (representing the "Read 2" reads of each read-pair).

Try <Execute command="ls reads_R1.fq.gz" inline="true" /> and <Execute command="zcat reads_R1.fq.gz | head -n 8" inline="true" /> to take a peek at the "R1" reads file (why did we pick a multiple of 4 in our `head` command?).

**2. The reference genome**

In the following steps, you'll map those reads to the <Link href="https://www.ncbi.nlm.nih.gov/nuccore/1798174254">SARS-CoV-2 reference genome</Link>.

We preloaded the reference genome's FASTA file, and its location is stored in the variable `$REF_FASTA`. Use <Execute command="echo $REF_FASTA" inline="true" /> to see the location.

You'll also need the reference genome's annotation GFF file, which is stored in the variable `$REF_GFF`. Use <Execute command="echo $REF_GFF" inline="true" /> to see the location.


**3. The primers**

For amplicon sequence data analysis, you'll also need a BED file representing the positions of the primers that were used in the amplicon sequencing protocol (we'll talk about these later in the tutorial), which is stored in the variable `$PRIMER_BED`. Use <Execute command="echo $PRIMER_BED" inline="true" /> to see its location.