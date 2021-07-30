<script>
import Execute from "../Execute.svelte";
import IGV from "../IGV.svelte";

let isOpen = false;
let igvOptions = {
	// Explicitly specify the reference URLs so that igv.js doesn't try downloading RefSeq genes
	reference: {
		id: "hg19",
		name: "Human (hg19)",
		fastaURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta",
		indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta.fai",
		cytobandURL: "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/cytoBand.txt",
		tracks: []
	},
	locus: "chr1:262,932-512,931",
	tracks: [
		{ url: "./data/cpg.bed", name: "CpG islands" },
		{ url: "./data/exons.bed", name: "RefSeq Exons" },
		{ url: "./data/gwas.bed", name: "GWAS SNPs" },
		{ url: "./data/hesc.chromHmm.bed", name: "chromHMM Predictions" }
	]
};
</script>

<p>We preloaded data from <a href="https://science.sciencemag.org/content/337/6099/1190" target="_blank">Maurano et al.</a> into your sandbox. Let's take a look at what files we now have:</p>

<p>Type <Execute command={"ls"} /> in the command line.</p>

<p>Your directory contains 7 <code>BED</code> files and 1 genome file. Three of these files (those starting with <code>f</code> for "fetal tissue") reflect Dnase I hypersensitivity sites measured in different fetal tissue samples.</p>

<p>In order to have a rough sense of the remaining <code>.bed</code> files, let's load them into IGV: <button class="btn btn-sm btn-primary" on:click={() => isOpen = !isOpen}>Launch IGV</button></p>

<IGV options={igvOptions} bind:isOpen={isOpen}>
	<span slot="after">
		Note that: 
		<ul>
			<li><code>cpg.bed</code> represents CpG islands in the human genome</li>
			<li><code>exons.bed</code> represents RefSeq exons from human genes</li>
			<li><code>gwas.bed</code> represents human disease-associated SNPs that were identified in genome-wide association studies (GWAS)</li>
			<li><code>hesc.chromHmm.bed</code> represents the predicted function (by chromHMM) of each interval in the genome of a human embryonic stem cell based upon ChIP-seq experiments from ENCODE</li>
		</ul>
	</span>	
</IGV>
