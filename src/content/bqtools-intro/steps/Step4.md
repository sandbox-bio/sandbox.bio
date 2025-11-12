<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

One of the benefits of BINSEQ is that we don't need to keep separate files for the different read pairs.
BINSEQ supports R1 and R2 pairing in the same file.

To create a paired file we can pass two files to the `bqtools encode` command:

<Execute command="bqtools encode fastq/sample1_R1.fastq.gz fastq/sample1_R2.fastq.gz"/>

By default, `bqtools` will autodetermine the output path from the sample names (in this case: `fastq/sample1.vbq`) - but you can override this with the `-o/--output` option.

<Execute command="bqtools encode fastq/sample1_R1.fastq.gz fastq/sample1_R2.fastq.gz -o my_output.vbq"/>
