<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

We will perform k-mer counting on two genomes: the Dengue virus and the Chikungunya virus, both transmitted by mosquitoes.

Their reference genomes has been preloaded into your environment:

<Execute command={`ls`} />

Let's start by counting k-mers in the Dengue genome using `k = 7` (we'll try other values later):

<Execute command={`jellyfish count \\ -m 7 \\ -s 11000 \\ -o dengue.jf \\ dengue.fa`} />

where:

- `-m 7`: count 7-mers
- `-s 11000`: our estimate of the number of k-mers in this genome. The Dengue genome is 10.6kbp so we choose 11k. This is just an approximation, Jellyfish will automatically resize its hash data structure to fit if needed
- `-o dengue.jf`: store results in `dengue.jf` (this is a binary file that we will query with Jellyfish)

<Alert>
   You can use `seqtk` to calculate how many bases are in a FASTA file:

<Execute command={`seqtk comp dengue.fa`} />

Use <Execute inline command="seqtk comp" /> to view what the columns in the output above represent.
</Alert>

Using the output of Jellyfish, let's look at the distribution of 7-mers:

<Execute command={`jellyfish histo dengue.jf`} />

This means that 4,264 k-mers are unique: they are only seen once in the Dengue genome. On the other hand, 1,538 k-mers are seen twice in the genome and 571 are seen three times.

<hr />

For k-mers that don't map to as many places in the genome, we can increase the value of `k`.

For example, let's count 9-mers:

<Execute command={`jellyfish count \\ -m 9 \\ -s 11000 \\ -o dengue.jf \\ dengue.fa`} />

Now you can see that the vast majority of 9-mers are unique within the Dengue genome.

<Execute command={`jellyfish histo dengue.jf`} />
