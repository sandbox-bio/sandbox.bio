<script>
import Execute from "$components/Execute.svelte";
</script>

With the same pile-up file from which called variants, we'll use `ivar` to call a "consensus sequence": a viral genome sequence containing the "consensus" (i.e., most abundant) nucleotide at each position. Run the following:

<Execute command="cat pileup.txt | \ ivar consensus \ -p consensus.fa \ -m 10 -t 0.5 -n N" />

Let's break this seemingly complex command into its individual components to make some sense of it:

- We're using the pipe character `|` to pipe the output of `cat` (via standard output) to `ivar` (via standard input)
- The `ivar` command we're piping the pile-up stream into is the following: `ivar consensus -p consensus.fa -m 10 -n N -t 0.5`
  - `consensus` tells `ivar` that we want to use its consensus-sequence-calling functionality
  - `-p consensus.fa` tells `ivar` that we want to write the output consensus genome sequence to the file `consensus.fa`
    - The output file is a FASTA file
  - `-m 10` tells `ivar` that we want to only call nucleotides in positions that have a **m**inimum depth of 10
    - In other words, positions with a coverage of less than 10 are ambiguous (i.e., we don't know what nucleotide is at that position)
  - `-t 0.5` tells `ivar` that we want to only call nucleotides in positions in which the most abundant nucleotide has a frequency of at least 0.5
    - In other words, positions in which the most abundant nucleotide appears in fewer than 50% of the reads covering that position are ambiguous
  - `-n N` tells `ivar` that we want to use `N` as the letter to represent ambiguous nucleotides
    - In other words, positions that are ambiguous (see `-m 10` and `-t 0.5` above) will contain `N` (rather than `A`, `C`, `G`, or `T`)

After running the above command, we will have successfully called a consensus genome sequence and written the resulting sequence to the file `consensus.fa`.

To see the contents of the consensus genome sequence output file, run the following:

<Execute command="cat consensus.fa" />

This consensus genome FASTA file is the other key result of our amplicon sequencing analysis, and it is what will be used in any downstream analyses (e.g. transmission clustering, phylogenetic inference, lineage assignment, etc.).
