<script>
import Execute from "$components/Execute.svelte";
</script>

Merging results in a new set of intervals representing the merged set of intervals in the input. That is, if a base pair in the genome is covered by 10 features, it will now only be represented once in the output file.

<Execute command={"bedtools merge -i exons.bed | \\ head -n 20"} />

A more sophisticated approach would be to not only merge overlapping intervals, but also report the **number** of intervals that were integrated into the new, merged interval. One does this with the `-c` and `-o` options. The `-c` option allows one to specify a column or columns in the input that you wish to summarize. The `-o` option defines the operation(s) that you wish to apply to each column listed for the `-c` option.  For example, to count the number of overlapping intervals that led to each of the new "merged" intervals, one will "count" the first column (though the second, third, fourth, etc. would work just fine as well).

<Execute command={"bedtools merge -i exons.bed -c 1 -o count | \\ head -n 20"} />

Many times you want to keep track of the details of exactly which intervals were merged. One way to do this is to create a list of the names of each feature. We can do with with the `collapse` operation available via the `-o` argument. The name of the exon is in the fourth column, so we ask `merge` to create a list of the exon names with `-c 4 -o collapse`:

<Execute command={"bedtools merge -i exons.bed -d 90 -c 1,4 -o count,collapse | \\ head -n 20"} />
