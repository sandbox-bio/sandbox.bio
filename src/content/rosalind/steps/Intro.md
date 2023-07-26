<script>
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
</script>

<Link href="http://rosalind.info"><img class="mb-3" src="/screenshot-rosalind.png" alt="Rosalind.info logo"></Link>

In this playground, you can write Python code to solve <Link href="http://rosalind.info">Rosalind.info</Link> exercises.

Go ahead, give it a try! Write some Python in the code editor on your right and press `Cmd + Enter` (or `Ctrl + Enter` on Windows) to run your code.

<Alert>
This playground uses <Link href="https://pyodide.org/">Pyodide</Link> to run a Python interpreter inside your browser!
</Alert>
