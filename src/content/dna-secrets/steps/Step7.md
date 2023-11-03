<script>
/*
	bowtie2 -x $REF -U reads.fq -S aligned.sam; samtools sort -o aligned.bam aligned.sam;  bcftools mpileup -f $REF_FASTA aligned.bam | bcftools call -m -v -Ob -o variants.bcf -; bcftools index variants.bcf

	bowtie2 -x $REF -U morereads.fq -S aligned2.sam; samtools sort -o aligned2.bam aligned2.sam;  bcftools mpileup -f $REF_FASTA aligned2.bam | bcftools call -m -v -Ob -o variants2.bcf -; bcftools index variants2.bcf

	bcftools merge variants.bcf variants2.bcf > combined.vcf
	
	bcftools query -f "%ALT" combined.vcf > secret
*/

import { onDestroy, onMount } from "svelte";
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
import Exercise from "$components/Exercise.svelte";
import { cli } from "$stores/cli";

// State
let dnaEncoded = "-";
let dnaDecoded = "";
let timers = [];
$: dnaDecoded = binaryToString(dnaEncoded.replaceAll("\n", "").split("").map(b => {
	// https://science.sciencemag.org/content/337/6102/1628
	if(b == "A" || b == "C") return "0";
	else return "1";
}).join("")) || "-";

// Converter
// https://stackoverflow.com/a/53247859
function binaryToString(input) {
	let bytesLeft = input;
	let result = '';

	// Check if we have some bytes left
	while (bytesLeft.length) {
		// Get the first digits
		const byte = bytesLeft.substr(0, 8);
		bytesLeft = bytesLeft.substr(8);
		result += String.fromCharCode(parseInt(byte, 2));
	}
	return result;
}

let criteria = [
{
	name: "File <code>secret</code> is a file containing the DNA secret",
	checks: [{
		type: "file",
		path: "secret",
		action: "contents",
		commandExpected: 'bcftools query -f "%ALT" combined.vcf'
	}]
}];

onMount(getSecret);
onDestroy(() => {
	timers.forEach((timer) => clearTimeout(timer));
});

async function getSecret() {
	try {
		const buffer = await $cli.readFile("/root/tutorial/secret");
		dnaEncoded = new TextDecoder().decode(buffer);
	} catch (error) {
		console.error(error);
	}
	timers.push(setTimeout(getSecret, 500));
}
</script>

Finally, it's time to decode the secret message!

<div class="form-floating mb-3">
	<input type="text" class="form-control" id="floatingInput" bind:value={dnaEncoded} disabled>
	<label for="floatingInput">DNA Sequence</label>
</div>
<div class="form-floating mb-3">
	<input type="text" class="form-control" id="floatingInput2" value={dnaDecoded} disabled>
	<label for="floatingInput2">Decoded Message</label>
</div>

Use the `bcftools query` command we introduced earlier to extract the `%ALT` column from the file `combined.vcf`, but show all the SNPs on one line (i.e. we don't need the `\n`).

<Exercise {criteria} />
