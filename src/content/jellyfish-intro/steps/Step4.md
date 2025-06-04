<script>
import Execute from "$components/Execute.svelte";
</script>

So far, we counted k-mers on entire genomes. In this section, we will count k-mers of sequencing reads and see if we can classify them to one or the other genome.

Let's simulate some sequencing reads from Dengue using `wgsim`:

<Execute command={`wgsim \\ -N 1000 \\ dengue.fa \\ r1.fastq \\ r2.fastq`} />

Now let's count k-mers in the sequencing reads, but only if they are also found in the Dengue genome (we'll ignore R2 reads):

<Execute command={`jellyfish count \\ -m 9 \\ -s 15000 \\ -C \\ -o map_to_dengue.jf \\ --if dengue.fa \\ r1.fastq`} />

Note that we're using `-C`, which means Jellyfish will only count **canonical k-mers**.

For example, the k-mers `ACCT` and `AGGT` are reverse complements of each other, so their counts are grouped under the so-called **canonical k-mer** `ACCT`, because it comes first alphabetically.

Note that we did not use `-C` earlier when counting k-mers on genomes. The difference is: when sequencing DNA, we can't differentiate which of the 2 strands the read comes from.

<hr />

Next, let's count k-mers from Dengue sequencing reads that are also found in the Chikungunya genome:

<Execute command={`jellyfish count \\ -m 9 \\ -s 15000 \\ -C \\ -o map_to_chikungunya.jf \\ --if chikungunya.fa \\ r1.fastq`} />

Looking at the distribution of k-mers, there are thousands more k-mers from Chikungunya that don't show up in the sequencing reads:

<Execute command={`jellyfish histo map_to_dengue.jf`} />

<Execute command={`jellyfish histo map_to_chikungunya.jf`} />

This is a simplistic example but illustrates how k-mers and their distribution can be a useful way to summarize sequencing reads for classification.
