<script>
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

To generate a simple QC report, we call `fastp` and specify the location of our paired-end data:

<Execute command="fastp \ --in1 ./HG004_R1.fastq.gz \ --in2 ./HG004_R2.fastq.gz" />

The structure of the output is as follows:

1. First, there are stats about how many reads and bases are in your dataset before and after filtering is applied. Note that the R1 dataset *decreases* from 25,000 reads to 24,455—i.e. `fastp` removed 545 reads that it deemed low quality). Also notice that the Q30 (number of bases with a quality > 30) *decreases* from 3,397,493 to 3,312,555—i.e. `fastp` removed 84,938 bases, and so the %Q30 *increases*.

2. Next, you'll see stats about how many reads and bases were removed and for what reason. For example, 942 reads were trimmed because adapter sequences were detected in those reads.

<Alert>Although the report mentions reads being filtered, `fastp` does not overwrite your original data! Later in this tutorial, we'll see how we can retrieve those filtered reads so we can use them in downstream analyses.</Alert>
