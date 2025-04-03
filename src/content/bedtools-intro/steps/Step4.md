<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

By default, `intersect` reports intervals that represent overlaps between your two files. To demonstrate, let's identify all of the CpG islands that overlap exons:

<Execute command={"bedtools intersect -a cpg.bed -b exons.bed | head -n 5"} />

<Alert>
	In this case, the intervals reported are NOT the original CpG intervals, but rather a refined interval reflecting solely the portion of each original CpG interval that overlapped with the exon(s).
</Alert>

The `-wa` (write A) and `-wb` (write B) options allow one to see the original records from the A and B files that overlapped.

As such, instead of not only showing you **where** the intersections occurred, it shows you **what** intersected.

<Execute command={"bedtools intersect -a cpg.bed -b exons.bed -wa -wb | head -n 5"} />
