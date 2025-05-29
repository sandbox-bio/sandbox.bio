<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

To run the split _k_-mer alignment, we first need to build an index of the split _k_-mers from the assemblies. This is accomplished with the `ska build` command which requires a tab-separated list consisting of the sequence name and assembly file path as input:

```bash
ID1	path/to/fasta_file_1.fasta.gz
ID2	path/to/fasta_file_2.fasta.gz
ID3	path/to/fasta_file_3.fasta.gz
.
.
.
```

which is provided in `assemblies/ska_input.tsv`. Let's first create an output folder with

<Execute command="mkdir output" />

Now, the split _k_-mers index is built from the input list by executing

<Execute command={"ska build \ -f assemblies/ska_input.tsv \ -k 31 \ -o output/ska_index"} />

which will take a bit of time to run. In a normal computer this would be very quick, but as we are running from a browser it will last 3-4 minutes or so.

In the example we chose the value of _k_ as 31 which is a commonly used choice for analysing bacterial genomes. For a more detailed analysis of the different possible values of _k_, have a look at the paper by <Link href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0258693">Bussi, Kapon, and Reich</Link> in PLOS ONE.
