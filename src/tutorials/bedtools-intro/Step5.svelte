<script>
import Execute from "../Execute.svelte";
</script>

<p>The <code>-wa</code> (write A) and <code>-wb</code> (write B) options allow one to see the original records from the A and B files that overlapped.</p>

<p>
	As such, instead of not only showing you <strong>where</strong> the intersections occurred, it shows you <strong>what</strong> intersected.
	<Execute command={"bedtools intersect -a cpg.bed -b exons.bed -wa -wb | head -n 5"} />
</p>
