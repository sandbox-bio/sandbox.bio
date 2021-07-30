<script>
import Execute from "../../Execute.svelte";
import Image from "../../Image.svelte";
</script>

<p>
	We often want to know which intervals of the genome are **NOT** "covered" by intervals in a given feature file. For example, if you have a set of ChIP-seq peaks, you may also want to know which regions of the genome are not bound by the factor you assayed. The `complement` addresses this task.
</p>

<Image src="https://bedtools.readthedocs.io/en/latest/_images/complement-glyph.png" alt="How bedtools complement works" />

<p></p>

<p>
	As an example, let's find all of the non-exonic (i.e., intronic or intergenic) regions of the genome.  Note, to do this you need a <a href="http://bedtools.readthedocs.org/en/latest/content/general-usage.html#genome-file-format" target="_blank">"genome" file</a>, which tells `bedtools` the length of each chromosome in your file. <i>Consider why the tool would need this information...</i>
</p>

<Execute command={"bedtools complement -i exons.bed -g genome.txt | head"} />
