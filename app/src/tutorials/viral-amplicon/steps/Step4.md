<script>
import Link from "components/Link.svelte";
import Execute from "components/Execute.svelte";
</script>

You may have noticed that, for the purposes of this exercise, the SAM file was just an intermediate file that was output by `minimap2` and then used as input by `viral_consensus`. It turns out that we can actually skip writing an intermediate file and instead directly feed the output of `minimap2` into `viral_consensus`. Run the following:

<Execute command="minimap2 -a -x sr $REF_FASTA \ reads_R1.fq reads_R2.fq | \ viral_consensus -i - \ -r $REF_FASTA \ -o consensus.fa \ -op position_counts.tsv \ -oi insertion_counts.json \ -p $PRIMER_BED -po 5" />

Let's break this seemingly complex command into its individual components to make some sense of it:

- First, we're calling the following command: `minimap2 -a -x sr $REF_FASTA reads_R1.fq reads_R2.fq`
  - This is essentially the same exact `minimap2` command as before
  - However, rather than having `minimap2` output the SAM results to a file (using the `-o FILE` argument), we're having it print the results to <Link href="https://en.wikipedia.org/wiki/Standard_streams#Standard_output_(stdout)">standard output</Link>
- Then, we're using the pipe character `|` to pipe the output of `minimap2` (via standard output) to `viral_consensus` (via <Link href="https://en.wikipedia.org/wiki/Standard_streams#Standard_input_(stdin)">standard input</Link>)
  - In general, the pipeline syntax `A | B` means "execute command `A`, and pipe whatever it outputs (via standard output) into command `B` (via standard input)"
  - Piping is useful for many reasons, e.g. minimizing slow disk access, reducing how many intermediate files we need to store (saving space), etc.
- The `viral_consensus` command we're piping the SAM stream into is the following: `viral_consensus -i - -r $REF_FASTA -o consensus.fa -op position_counts.tsv -oi insertion_counts.json -p $PRIMER_BED -po 5`
  - This is essentially the same exact `viral_consensus` command as before
  - However, rather than having `viral_consensus` read an input file (using the `-i FILE` argument), we're having it read the input data from standard input (by specifying `-i -`)

After running the above command, we will have performed the exact same analysis as before, but with all data processing happening on-the-fly, without the need to save an intermediate SAM file.
