<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

Let's output the first few k-mers and their counts in FASTA format:

<Execute command={`jellyfish dump dengue.jf | head`} />

For example, this FASTA record:

```
  >2
  AAGTTTTCA
```

means the k-mer `AAGTTTTCA` was seen twice.

<hr />

To query for a particular k-mer of interest, say `ACAGTGGAC`, you can use `jellyfish query`:

<Execute command={`jellyfish query dengue.jf ACAGTGGAC`} />

<Execute command={`jellyfish query chikungunya.jf ACAGTGGAC`} />

This tells us that the `ACAGTGGAC` k-mer is found in Denge but not Chikungunya.

<hr />

To get k-mers found in Dengue but not Chikungunya, we can use `jellyfish count --if`:

<Execute command={`jellyfish count \\ -m 9 \\ -s 15000 \\ -o intersect.jf \\ --if chikungunya.fa \\ dengue.fa`} />

The distribution of k-mers now looks different:

<Execute command={`jellyfish histo intersect.jf`} />

This means there are 10,764 k-mers from the Chikungunya genome that were _not_ found in the Dengue genome.
