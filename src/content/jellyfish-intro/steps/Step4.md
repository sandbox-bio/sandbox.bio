<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

So far, we counted k-mers on entire genomes. In this section, we will count k-mers of sequencing reads and see if we can classify them to one or the other genome.

Let's simulate some sequencing reads from Dengue using `wgsim`:

<Execute command={`wgsim -N 1000 dengue.fa r1.fastq r2.fastq`} />

We can count k-mers from the sequencing reads only if they are found in the Dengue genome (we'll ignore R2 reads here):

<Execute command={`jellyfish count -m 9 -s 15000 -C -o map_to_dengue.jf --if dengue.fa r1.fastq`} />

Note that we're using `-C`, which means Jellyfish will only count **canonical k-mers**. For example, the 4-mers `ACCT` and `AGGT` are reverse complements of each other. Their counts will be grouped under the so-called k-mer `ACCT`, which is the canonical one since it comes first in alphabetically. Note that we did not use `-C` earlier when counting k-mers on genomes. The difference is: when sequencing DNA, we can't differentiate which of the 2 strands the read comes from.

Repeat the same command, but for Chikungunya:

<Execute command={`jellyfish count -m 9 -s 15000 -C -o map_to_chikungunya.jf --if chikungunya.fa r1.fastq`} />

Now, looking at the distribution of k-mers, there are thousands more k-mers from Chikungunya that don't show up in the sequencing reads:

<Execute command={`jellyfish histo map_to_dengue.jf`} />

<Execute command={`jellyfish histo map_to_chikungunya.jf`} />

This is a simplistic example but illustrates how k-mers and their distribution can be a useful way to summarize sequencing reads for classification.
