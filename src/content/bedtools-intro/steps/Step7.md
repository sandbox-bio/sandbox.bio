<script>
import Execute from "$components/Execute.svelte";
</script>

As of version `2.21.0`, bedtools is able to intersect an "A" file against one or more "B" files. This greatly simplifies analyses involving multiple datasets relevant to a given experiment. For example, let's intersect exons with CpG islands, GWAS SNPs, an the ChromHMM annotations:

<Execute command="bedtools intersect -a exons.bed -b cpg.bed gwas.bed hesc.chromHmm.bed -sorted | head" />

Now by default, this isn't incredibly informative as we can't tell which of the three "B" files yielded the intersection with each exon. However, if we use the `-wa` and `wb` options, we can see from which file number (following the order of the files given on the command line) the intersection came. In this case, the 7th column reflects this file number:

<Execute command="bedtools intersect -a exons.bed -b cpg.bed gwas.bed hesc.chromHmm.bed -sorted -wa -wb \\ | head -n 10000 \\ | tail -n 10" />

Additionally, one can use file "labels" instead of file numbers to facilitate interpretation, especially when there are _many_ files involved:

<Execute command="bedtools intersect -a exons.bed -b cpg.bed gwas.bed hesc.chromHmm.bed -sorted -wa -wb -names cpg gwas chromhmm \\ | head -n 10000 \\ | tail -n 10" />
