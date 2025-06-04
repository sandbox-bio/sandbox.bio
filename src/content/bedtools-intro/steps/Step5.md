<script>
import Execute from "$components/Execute.svelte";
</script>

The `-wo` (write overlap) option allows one to also report the **number** of base pairs of overlap between the features that overlap between each of the files:

<Execute command="bedtools intersect -a cpg.bed -b exons.bed -wo | head" />

We can also count, for each feature in the "A" file, the number of overlapping features in the "B" file. This is handled with the `-c` option:

<Execute command="bedtools intersect -a cpg.bed -b exons.bed -c | head" />

Often we want to identify those features in our A file that **do not** overlap features in the B file. The `-v` option is your friend in this case:

<Execute command="bedtools intersect -a cpg.bed -b exons.bed -v | head" />

Recall that the default is to report overlaps between features in A and B so long as **at least one basepair** of overlap exists. However, the `-f` option allows you to specify what fraction of each feature in A should be overlapped by a feature in B before it is reported. Let's be more strict and require 50% of overlap:

<Execute command="bedtools intersect -a cpg.bed -b exons.bed -wo -f 0.50 | head" />
