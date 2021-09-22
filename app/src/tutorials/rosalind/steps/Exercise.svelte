<script>
import { Converter } from "showdown";
import Alert from "./components/Alert.svelte";
import { tutorial } from "./stores/tutorials";

const converter = new Converter();

$: console.log($tutorial.step);

let problem = {"id": "RNA", "title": "Transcribing DNA into RNA", "url": "http://rosalind.info/problems/rna/", "given": "A DNA string _t_ having length at most 1000 nt.", "return": "The transcribed RNA string of _t_.", "sample_data": "GATGGAACTTGACTACGTAAATT\\n", "sample_output": "GAUGGAACUUGACUACGUAAAUU\\n"};
</script>

<Alert>
	Submit your answer on <a href={problem.url} target="_blank">Rosalind</a>
</Alert>

<h6>Given:</h6>
{@html converter.makeHtml(problem.given)}

<h6>Return:</h6>
{@html converter.makeHtml(problem.return)}

<h6>Sample Input:</h6>
{@html converter.makeHtml(problem.sample_data)}

<h6>Sample Output:</h6>
{@html converter.makeHtml(problem.sample_output)}

<!-- ##### Problem

**Given**: {problem.given}

**Return**: The transcribed RNA string of _t_.

##### Sample Dataset

```
GATGGAACTTGACTACGTAAATT
```

##### Sample Output

```
GAUGGAACUUGACUACGUAAAUU
``` -->
