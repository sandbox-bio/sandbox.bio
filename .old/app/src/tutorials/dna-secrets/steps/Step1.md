<script>
import Link from "components/Link.svelte";
import Execute from "components/Execute.svelte";
</script>

**1. The sequencing reads**

The reads have been preloaded to your working directory. Try <Execute command="ls reads.fq" inline /> and <Execute command="head -n 8 reads.fq" inline /> (why did we pick a multiple of 4 in our `head` command?).

**2. The reference genome**

In the following steps, you'll map those reads to the <Link href="https://en.wikipedia.org/wiki/Lambda_phage">Lambda phage</Link> reference genome.

We preloaded the `bowtie2` indexes for that reference genome; its location is stored in the variable `$REF`. Use <Execute command="echo $REF" inline /> to see the location.

For variant calling, you'll also need the reference genome's FASTA file, and its location is stored in the variable `$REF_FASTA`. Use <Execute command="echo $REF_FASTA" inline /> to see the location (or use <Execute command="env" inline /> to list all variables).
