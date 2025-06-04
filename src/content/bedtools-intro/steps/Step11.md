<script>
import Execute from "$components/Execute.svelte";
import Image from "$components/Image.svelte";
import Link from "$components/Link.svelte";
</script>

We often want to know which intervals of the genome are **NOT** "covered" by intervals in a given feature file. For example, if you have a set of ChIP-seq peaks, you may also want to know which regions of the genome are not bound by the factor you assayed. The `complement` addresses this task.

<Image src="/data/bedtools-intro/complement-glyph.png" alt="How bedtools complement works" />

As an example, let's find all of the non-exonic (i.e., intronic or intergenic) regions of the genome. Note, to do this you need a <Link href="http://bedtools.readthedocs.org/en/latest/content/general-usage.html#genome-file-format">genome file</Link>, which tells `bedtools` the length of each chromosome in your file. _Consider why the tool would need this information..._

<Execute command="bedtools complement -i exons.bed -g genome.txt | head" />
