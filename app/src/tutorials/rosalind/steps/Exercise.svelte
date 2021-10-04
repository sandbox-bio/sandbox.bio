<script>
import { Converter } from "showdown";
import Alert from "./components/Alert.svelte";
import { tutorial } from "./stores/tutorial";

const converter = new Converter();

$: exercise = $tutorial.steps[$tutorial.step].rosalind;
</script>

<Alert>
	Submit your answer on <a href="http://rosalind.info/problems/{exercise.id.toLowerCase()}" target="_blank">Rosalind</a>
</Alert>

<h6>Given:</h6>
{@html converter.makeHtml(exercise.given.trim())}

<h6>Return:</h6>
{@html converter.makeHtml(exercise.return.trim())}

<h6>Sample Input:</h6>
<code class="text-break">
	{@html exercise.sample_data.replaceAll("\n", "<br>")}
</code>

<h6 class="mt-3">Sample Output:</h6>
<code class="text-break">
	{@html exercise.sample_output.replaceAll("\n", "<br>")}
</code>

<style>
code {
	font-size: 0.8em;
	color: #d63384;
}
</style>
