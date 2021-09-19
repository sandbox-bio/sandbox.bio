<script>
import Alert from "./components/Alert.svelte";
import { tutorial } from "./stores/tutorials";

$: console.log($tutorial.step);
</script>

<Alert>
	Submit your answer on <a href="#" target="_blank">Rosalind</a>
</Alert>


