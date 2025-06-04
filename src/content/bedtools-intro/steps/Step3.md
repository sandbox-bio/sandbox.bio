<script>
import Image from "$components/Image.svelte";
</script>

The `intersect` command is the workhorse of the `bedtools` suite. It compares two or more BED/BAM/VCF/GFF files and identifies all the regions in the gemome where the features in the two files overlap (that is, share at least one base pair in common).

<Image src="/data/bedtools-intro/intersect-glyph.png" alt="How bedtools intersect works with one or more files" />
