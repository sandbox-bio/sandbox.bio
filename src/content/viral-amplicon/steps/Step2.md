<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

First, we'll use `minimap2` to align the sequencing paired-end reads in `reads_R1.fq` and `reads_R2.fq`to the reference genome located at `$REF`. While `minimap2` aligns the reads, we will <Link href="https://en.wikipedia.org/wiki/Pipeline_(Unix)#Pipelines_in_command_line_interfaces">pipe</Link> the resulting SAM stream to `samtools` to sort. Run the following:

<Execute command="minimap2 -a -o reads.mapped.sam -x sr \ $REF_FASTA reads_R1.fq reads_R2.fq" />

Let's make some sense of this `minimap2` command:
- `-a` tells `minimap2` that we want it to output the results in the SAM format
- `-o reads.mapped.sam` tells `minimap2` that we want to write the output SAM results to the file `reads.mapped.sam` (rather than printing them to standard output)
- `-x sr` tells `minimap2` that we want to use its **s**hort-**r**ead mapping preset
  - This is because our reads were sequenced using an Illumina short-read sequencer
- Next, we're specifying the reference genome: `$REF_FASTA`
- Lastly, we're specifying the read files we want to map: `reads_R1.fq` and `reads_R2.fq`

After running the above command, we will have successfully mapped our reads (`reads_R1.fq` and `reads_R2.fq`) to the reference genome (`$REF`) and written the results to the file `reads.mapped.bam`.

To see the first few lines of the SAM output file in the human-readable SAM format, run the following:

<Execute command="head -n 5 reads.mapped.sam" />

The first few lines (beginning with `@`) are SAM header lines, and the rest of the lines are SAM alignments, one line per read or mate. See the <Link href="http://samtools.sourceforge.net/SAM1.pdf">SAM specification</Link> for details about how to interpret the SAM file format.
