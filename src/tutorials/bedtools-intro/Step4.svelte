<script>
import Execute from "../Execute.svelte";
</script>

<p>
	By default, <code>intersect</code> reports intervals that represent overlaps between your two files. To demonstrate, let's identify all of the CpG islands that overlap exons:

	<Execute command={"bedtools intersect -a cpg.bed -b exons.bed | head -n 5"} />
</p>

<div class="alert alert-info m-2">
	<strong>Note:</strong> In this case, the intervals reported are NOT the original CpG intervals, but rather a refined interval reflecting solely the portion of each original CpG interval that overlapped with the exon(s).
</div>
