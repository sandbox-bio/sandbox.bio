<script>
import Execute from "$components/Execute.svelte";
import IGVModal from "$components/igv/IGVModal.svelte";
import Link from "$components/Link.svelte";

let isOpen = false;
let igvOptions = {
	locus: "chr1:262,932-512,931",
	tracks: [
		{ url: "/data/bedtools-intro/cpg.bed", name: "CpG islands" },
		{ url: "/data/bedtools-intro/exons.bed", name: "RefSeq Exons" },
		{ url: "/data/bedtools-intro/gwas.bed", name: "GWAS SNPs" },
		{ url: "/data/bedtools-intro/hesc.chromHmm.bed", name: "chromHMM Predictions" }
	]
};
</script>

We preloaded data from <Link href="https://science.sciencemag.org/content/337/6099/1190">Maurano et al.</Link> into your sandbox. Let's take a look at what files we now have:

Type <Execute command="ls" inline /> in the command line.

Your directory contains 7 `BED` files and 1 genome file. Three of these files (those starting with `f` for "fetal tissue") reflect Dnase I hypersensitivity sites measured in different fetal tissue samples.

In order to have a rough sense of the remaining `.bed` files, let's load them into IGV: <button class="btn btn-sm btn-primary" on:click={() => isOpen = !isOpen}>Launch IGV</button>

<IGVModal options={igvOptions} bind:isOpen={isOpen}>
	<span slot="after">
		Note that:

    	* `cpg.bed` represents CpG islands in the human genome
    	* `exons.bed` represents RefSeq exons from human genes
    	* `gwas.bed` represents human disease-associated SNPs that were identified in genome-wide association studies (GWAS)
    	* `hesc.chromHmm.bed` represents the predicted function (by chromHMM) of each interval in the genome of a human embryonic stem cell based upon ChIP-seq experiments from ENCODE
    </span>

</IGVModal>
