<script>
// Shortcut: samtools view -b sample.sam -o sample.bam; samtools sort sample.bam -o sample.sorted.bam; samtools index sample.sorted.bam; 

import Execute from "$components/Execute.svelte";
</script>

As we discussed earlier, the FLAG field in the BAM format encodes several key
pieces of information regarding how an alignment aligned to the reference genome.
We can exploit this information to isolate specific types of alignments that we
want to use in our analysis.

For example, we often want to call variants solely from paired-end sequences
that aligned "properly" to the reference genome.

To ask the `view` command to report solely "proper pairs", we use the `-f` option
and ask for alignments where the second bit is true (proper pair is true):

<Execute command={"samtools view -f 0x2 sample.sorted.bam | head"} />

How many _properly_ paired alignments are there? (use the `-c` option)

<Execute command={"samtools view -c -f 0x2 sample.sorted.bam"} />

Now, let's ask for alignments that are NOT properly paired. To do this, we use the `-F` option (note the capitalization to denote "opposite").

<Execute command={"samtools view -c -F 0x2 sample.sorted.bam"} />

How many _total_ alignments?

<Execute command={"samtools view -c sample.sorted.bam"} />

Does everything add up?

To get a summary of the flags in our BAM file, we can use `samtools flagstats`:

<Execute command={"samtools flagstats sample.sorted.bam"} />
