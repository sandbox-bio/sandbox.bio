<script>
import Execute from "../Execute.svelte";
</script>

<p>We preloaded data from <a href="https://science.sciencemag.org/content/337/6099/1190" target="_blank">Maurano et al.</a> into your sandbox. Let's take a look at what files we now have:</p>

<p>Type <Execute command={"ls"} /> in the command line.</p>

<p>Your directory contains 7 <code>BED</code> files and 1 genome file. Use <code>head</code> to preview these files: <Execute command={"head cpg.bed"} /></p>
	
<p>Three of these files (those starting with <code>f</code> for "fetal tissue") reflect Dnase I hypersensitivity sites measured in different fetal tissue samples.</p>

<span>
	In addition: 
	<ul>
		<li><code>cpg.bed</code> represents CpG islands in the human genome</li>
		<li><code>exons.bed</code> represents RefSeq exons from human genes</li>
		<li><code>gwas.bed</code> represents human disease-associated SNPs that were identified in genome-wide association studies (GWAS)</li>
		<li><code>hesc.chromHmm.bed</code> represents the predicted function (by chromHMM) of each interval in the genome of a human embryonic stem cell based upon ChIP-seq experiments from ENCODE</li>
	</ul>
</span>

<p>
The latter 4 files were extracted from the UCSC Genome Browser's <a href="http://genome.ucsc.edu/cgi-bin/hgTables?command=start" target="_blank">Table Browser</a>.
</p>

<p>
In order to have a rough sense of the data, let's load the `cpg.bed`, `exons.bed`, `gwas.bed`, and `hesc.chromHmm.bed` files into [IGV](http://www.broadinstitute.org/igv/). To do this, launch IGV, then click File->Load from File. Then select the four files. IGV will warn you that you need to create an index for a couple of the files. Just click OK, as these indices are created automatically and speed up the processing for IGV.
</p>

<p>
Visualization in IGV or other browsers such as UCSC is a tremendously useful way to make sure that your results make sense to your eye. Conveniently, a subset of bedtools is built-into IGV!
</p>
