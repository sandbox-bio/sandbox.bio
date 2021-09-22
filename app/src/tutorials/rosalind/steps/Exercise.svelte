<script>
import { Converter } from "showdown";
import Alert from "./components/Alert.svelte";
import { tutorial } from "./stores/tutorial";

$: exercise = $tutorial.steps[$tutorial.step].rosalind;

const converter = new Converter();
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
{@html converter.makeHtml(exercise.sample_data.trim())}
</code>

<h6>Sample Output:</h6>
<code class="text-break">
{@html converter.makeHtml(exercise.sample_output.trim())}
</code>
