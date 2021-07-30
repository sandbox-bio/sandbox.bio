<script>
import Execute from "../Execute.svelte";
import Image from "../Image.svelte";
</script>

<p>
	So far the examples presented have used the traditional algorithm in bedtools for finding intersections.  It turns out, however, that bedtools is much faster when using presorted data.
</p>

<p>
	For example, compare the difference in speed between the two approaches when finding intersections between <code>exons.bed</code> and <code>hesc.chromHmm.bed</code>:
</p>

<Execute command={"time bedtools intersect -a gwas.bed -b hesc.chromHmm.bed > /dev/null"} />

<p></p>

<Execute command={"time bedtools intersect -a gwas.bed -b hesc.chromHmm.bed -sorted > /dev/null"} />

<div class="p-2 mt-3">
	<div class="alert alert-info p-3">
		<strong>Note:</strong> While the run times in this example are quite small, the performance gains from using the <code>-sorted</code> option grow as datasets grow larger. For example, compare the runtimes of the sorted and unsorted approaches as a function of dataset size in the figure below. The important thing to
		remember is that each dataset must be sorted by chromosome and then by start position: <code>sort -k1,1 -k2,2n</code>.
	</div>
</div>

<Image alt="Performance gains due to using sorted data increase as the number of BAM alignments increases" src="https://bedtools.readthedocs.io/en/latest/_images/speed-comparo.png" />
