<script>
import { onMount } from "svelte";
import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";

// State
let dnaEncoded = "-";
let dnaDecoded = "";
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

onMount(async () => {
	setInterval(async () => {
		dnaEncoded = await $CLI.exec("cat secret");
	}, 700);
});
</script>

<Alert>
	**Note**: This problem set assumes you already completed the tutorials on <Link href="/tutorials?id=bowtie2-intro">bowtie2</Link> and <Link href="/tutorials?id=samtools-intro">samtools</Link>, but googling works too 🙂
</Alert>

In this problem set, you're given DNA sequencing reads that encode a secret message. To _decode_ this message, you'll have to align these reads to a reference genome and call variants. The SNPs you identify via variant calling, when sorted by genomic coordinate, will make up a string of DNA letters (`A`, `C`, `G`, `T`) that can be decoded to reveal the secret message.

The final goal of this tutorial is to create a file called `secret` that contains that DNA sequence and use our decoder to see what the underlying message is.

Here's a simple DNA decoder below based on the algorithm described in <Link href="https://science.sciencemag.org/content/337/6102/1628">Church et al, 2013</Link>. It will show the decoded value of the DNA stored in the file `secret`.

<div class="form-floating mb-3">
	<input type="text" class="form-control" id="floatingInput" bind:value={dnaEncoded} disabled>
	<label for="floatingInput">DNA Sequence</label>
</div>
<div class="form-floating mb-3">
	<input type="text" class="form-control" id="floatingInput2" value={dnaDecoded} disabled>
	<label for="floatingInput2">Decoded Message</label>
</div>

For example, try:

<Execute command='echo "CGGCGAACAGGCCTAGATTAGGCCCTTCTTCCCGGCGGTG" > secret' />
