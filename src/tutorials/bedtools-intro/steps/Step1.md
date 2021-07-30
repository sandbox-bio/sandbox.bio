<script>
import Execute from "../../Execute.svelte";
import IGV from "../../IGV.svelte";
import Link from "../../Link.svelte";

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

We preloaded data from <Link href="https://science.sciencemag.org/content/337/6099/1190">Maurano et al.</Link> into your sandbox. Let's take a look at what files we now have:

Type <Execute command={"ls"} /> in the command line.

Your directory contains 7 `BED` files and 1 genome file. Three of these files (those starting with `f` for "fetal tissue") reflect Dnase I hypersensitivity sites measured in different fetal tissue samples.

In order to have a rough sense of the remaining `.bed` files, let's load them into IGV: <button class="btn btn-sm btn-primary" on:click={() => isOpen = !isOpen}>Launch IGV</button>

<IGV options={igvOptions} bind:isOpen={isOpen}>
	<span slot="after">
		Note that: 

		* `cpg.bed` represents CpG islands in the human genome
		* `exons.bed` represents RefSeq exons from human genes
		* `gwas.bed` represents human disease-associated SNPs that were identified in genome-wide association studies (GWAS)
		* `hesc.chromHmm.bed` represents the predicted function (by chromHMM) of each interval in the genome of a human embryonic stem cell based upon ChIP-seq experiments from ENCODE
	</span>	
</IGV>
