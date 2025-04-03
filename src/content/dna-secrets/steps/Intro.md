<script>
import { onDestroy, onMount } from "svelte";
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
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

<Alert>
	Here we assume you completed the <Link href="/tutorials/bowtie2-intro">bowtie2</Link> and <Link href="/tutorials/samtools-intro">samtools</Link> tutorials, but googling works too 🙂
</Alert>

In this problem set, you're given DNA sequencing data that, when analyzed by mapping reads and calling variants, reveals a secret DNA message. Here's a decoder based on <Link href="https://science.sciencemag.org/content/337/6102/1628">Church et al.</Link> that decodes DNA sequences stored in the file `secret`:

<div class="form-floating mb-1">
	<input type="text" class="form-control" id="floatingInput" bind:value={dnaEncoded} disabled>
	<label for="floatingInput">DNA Sequence</label>
</div>
<div class="form-floating mb-3">
	<input type="text" class="form-control" id="floatingInput2" value={dnaDecoded} disabled>
	<label for="floatingInput2">Decoded Message</label>
</div>

For example, try:

- <Execute command='echo "ATGAGACTCTGGACTTCGTCGGGAAATAAGGTATGTATAA" > secret' />
- <Execute command='echo "CTGGCTCCCTTCGAAAAGGATCATAGTTAAGT" > secret' />
- <Execute command='echo "CTGAATTCCGGTCGCTCGGCTGTCACTGGTTG" > secret' />
