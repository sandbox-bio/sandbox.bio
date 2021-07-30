<script>
import Execute from "../Execute.svelte";
</script>

<p>
	The <code>-wo</code> (write overlap) option allows one to also report the <strong>number</strong> of base pairs of overlap between the features that overlap between each of the files:
</p>

<Execute command={"bedtools intersect -a cpg.bed -b exons.bed -wo | head"} />

<p></p><p></p>

<p>
	We can also count, for each feature in the "A" file, the number of overlapping features in the "B" file. This is handled with the <code>-c</code> option:
</p>

<Execute command={"bedtools intersect -a cpg.bed -b exons.bed -c | head"} />

<p></p><p></p>

<p>
	Often we want to identify those features in our A file that <strong>do not</strong> overlap features in the B file. The <code>-v</code> option is your friend in this case:
</p>

<Execute command={"bedtools intersect -a cpg.bed -b exons.bed -v | head"} />

<p></p><p></p>

<p>
	Recall that the default is to report overlaps between features in A and B so long as <strong>at least one basepair</strong> of overlap exists. However, the <code>-f</code> option allows you to specify what fraction of each feature in A should be overlapped by a feature in B before it is reported. Let's be more strict and require 50% of overlap:
</p>

<Execute command={"bedtools intersect -a cpg.bed -b exons.bed -wo -f 0.50 | head"} />
