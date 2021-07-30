<script>
import Execute from "../Execute.svelte";
</script>

<p>
The <code>intersect</code> command is the workhorse of the <code>bedtools</code> suite. It compares two or more BED/BAM/VCF/GFF files and identifies all the regions in the gemome where the features in the two files overlap (that is, share at least one base pair in common).
</p>

<img alt="How bedtools intersect works with one or more files" src="https://bedtools.readthedocs.io/en/latest/_images/intersect-glyph.png">

<p>
	By default, <code>intersect</code> reports intervals that represent overlaps between your two files. To demonstrate, let's identify all of the CpG islands that overlap exons:

	<Execute command={"bedtools intersect -a cpg.bed -b exons.bed | head -n 5"} />
</p>

<div class="alert alert-info m-2">
	<strong>Note:</strong> In this case, the intervals reported are NOT the original CpG intervals, but rather a refined interval reflecting solely the portion of each original CpG interval that overlapped with the exon(s).
</div>
