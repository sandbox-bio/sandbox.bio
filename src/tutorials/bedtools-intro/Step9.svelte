<script>
import Execute from "../Execute.svelte";
</script>

<p>Merging results in a new set of intervals representing the merged set of intervals in the input. That is, if a base pair in the genome is covered by 10 features, it will now only be represented once in the output file.</p>

<Execute command={"bedtools merge -i exons.bed | head -n 20"} />

<p></p><p></p>

<p>A more sophisticated approach would be to not only merge overlapping intervals, but also report the <strong>number</strong> of intervals that were integrated into the new, merged interval. One does this with the <code>-c</code> and <code>-o</code> options. The <code>-c</code> option allows one to specify a column or columns in the input that you wish to summarize. The <code>-o</code> option defines the operation(s) that you wish to apply to each column listed for the <code>-c</code> option.  For example, to count the number of overlapping intervals that led to each of the new "merged" intervals, one will "count" the first column (though the second, third, fourth, etc. would work just fine as well).</p>

<Execute command={"bedtools merge -i exons.bed -c 1 -o count | head -n 20"} />

<p></p><p></p>

<p>
	Many times you want to keep track of the details of exactly which intervals were merged. One way to do this is to create a list of the names of each feature. We can do with with the <code>collapse</code> operation available via the <code>-o</code> argument. The name of the exon is in the fourth column, so we ask <code>merge</code> to create a list of the exon names with <code>-c 4 -o collapse</code>:
</p>

<Execute command={"bedtools merge -i exons.bed -d 90 -c 1,4 -o count,collapse | head -n 20"} />
