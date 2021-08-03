<script>
import Execute from "components/Execute.svelte";
</script>

The samtools `view` command is the most versatile tool in the samtools package.
It's main function, not surprisingly, is to allow you to convert the binary
(i.e., easy for the computer to read and process) alignments in the BAM file
view to text-based SAM alignments that are easy for *humans* to read and process.

**Scrutinize some alignments**

Let us start by inspecting the first five alignments in our BAM in detail:

<Execute command={"samtools view sample.sorted.bam | head -n 5"} />

**Let's make the FLAG more readable**

You can use the detailed help to get a better sense of what each character in the human "readable" FLAG means:

<Execute command={"samtools flags"} />

**Count the alignments**

To count the total number of alignments:

<Execute command={"samtools view -c sample.sorted.bam"} />

**Inspect the header**

The "header" in a BAM file records important information regarding the 
reference genome to which the reads were aligned, as well as other information
about how the BAM has been processed. One can ask the `view`
command to report solely the header by using the `-H` option.

<Execute command={"samtools view -H sample.sorted.bam"} />

**Capture the FLAG**

As we discussed earlier, the FLAG field in the BAM format encodes several key
pieces of information regarding how an alignment aligned to the reference genome.
We can exploit this information to isolate specific types of alignments that we
want to use in our analysis.

For example, we often want to call variants solely from paired-end sequences
that aligned "properly" to the reference genome.

To ask the `view` command to report solely "proper pairs", we use the `-f` option
and ask for alignments where the second bit is true (proper pair is true):

<Execute command={"samtools view -f 0x2 sample.sorted.bam"} />

How many *properly* paired alignments are there?

<Execute command={"samtools view -c -f 0x2 sample.sorted.bam"} />

Now, let's ask for alignments that are NOT properly paired.  To do this, we use the `-F` option (note the capitalization to denote "opposite").

<Execute command={"samtools view -F 0x2 sample.sorted.bam"} />

How many *improperly* paired alignments are there?

<Execute command={"samtools -c view -F 0x2 sample.sorted.bam"} />

Does everything add up?
