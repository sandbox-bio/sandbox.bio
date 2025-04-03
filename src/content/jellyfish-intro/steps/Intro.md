<script>
import Link from "$components/Link.svelte";
</script>

K-mer counting is a tool used in many areas of bioinformatics, including metagenomics and genome assembly.

A k-mer is simply a substring of length `k`. For example, the sequence `ACGTA` has four 2-mers (k = 2):

```
ACGTA
AC
 CG
  GT
   TA
```

In this tutorial, we'll use the tool <Link href="https://github.com/gmarcais/Jellyfish">Jellyfish</Link> to count k-mers in genomes and sequencing reads.
