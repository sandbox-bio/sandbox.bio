<script>
import Link from "components/Link.svelte";
import Execute from "components/Execute.svelte";
</script>

First, we'll use `minimap2` to align the sequencing paired-end reads in `reads_R1.fq.gz` and `reads_R2.fq.gz`to the reference genome located at `$REF`. While `minimap2` aligns the reads, we will <Link href="https://en.wikipedia.org/wiki/Pipeline_(Unix)#Pipelines_in_command_line_interfaces">pipe</Link> the resulting SAM stream to `samtools` to sort. Run the following:

<Execute command="minimap2 -a -x sr $REF_FASTA \ reads_R1.fq.gz reads_R2.fq.gz | \ samtools sort -o untrimmed.sorted.bam" />

Let's break this seemingly complex command into its individual components to make some sense of it:

- First, we're calling the following command: `minimap2 -a -x sr $REF_FASTA reads_R1.fq.gz reads_R2.fq.gz`
  - `-a` tells `minimap2` that we want it to output the results in the SAM format
  - `-x sr` tells `minimap2` that we want to use its **s**hort-**r**ead mapping preset
    - This is because our reads were sequenced using an Illumina short-read sequencer
  - Next, we're specifying the reference genome: `$REF_FASTA`
  - Lastly, we're specifying the read files we want to map: `reads_R1.fq.gz` and `reads_R2.fq.gz`
  - Rather than having `minimap2` output the SAM results to a file (using the `-o FILE` argument), we're having it print the results to <Link href="https://en.wikipedia.org/wiki/Standard_streams#Standard_output_(stdout)">standard output</Link>
- Then, we're using the pipe character `|` to pipe the output of `minimap2` (via standard output) to `samtools` (via <Link href="https://en.wikipedia.org/wiki/Standard_streams#Standard_input_(stdin)">standard input</Link>)
  - In general, the pipeline syntax `A | B` means "execute command `A`, and pipe whatever it outputs (via standard output) into command `B` (via standard input)"
  - Piping is useful for many reasons, e.g. minimizing slow disk access, reducing how many intermediate files we need to store (saving space), etc.
- The `samtools` command we're piping the SAM stream into is the following: `samtools sort -o untrimmed.sorted.bam`
  - `sort` tells `samtools` that we want to use its sorting functionality
  - `-o untrimmed.sorted.bam` tells `samtools` that we want to write the output BAM results to the file `untrimmed.sorted.bam` (rather than printing them to standard output)
    - BAM is a binary (i.e., non-human-readable) and (potentially) compressed file format that stores the same information as the SAM format
  - Rather than having `samtools` read an input file (by specifying a positional argument at the end of the command), we're having it read the input data from standard input

After running the above command, we will have sucessfully mapped our reads (`reads_R1.fq.gz` and `reads_R2.fq.gz`) to the reference genome (`$REF`), sorted the results, and written the sorted results to the file `untrimmed.sorted.bam`.

To see the first few lines of the BAM output file in the human-readable SAM format, run the following:

<Execute command="samtools view -h untrimmed.sorted.bam | \ head -n 5" />

The first few lines (beginning with `@`) are SAM header lines, and the rest of the lines are SAM alignments, one line per read or mate. See the <Link href="http://samtools.sourceforge.net/SAM1.pdf">SAM specification</Link> for details about how to interpret the SAM file format.
