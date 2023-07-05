<script>
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

<Alert>
	Make sure you're comfortable with the <Link href="/tutorials?id=terminal-basics">Terminal Basics</Link> (e.g. `ls`, `head`, `tail`, `grep`) before going through this tutorial.
</Alert>

Awk is a power tool to help you **filter, extract and transform data files** on the command line.

It is most commonly used with tab- and comma-separated files, where the idea is to apply a certain transformation to every line (and sometimes column) of a file. As you'll see later in this tutorial, `awk` is actually a full-fledged programming language!

Here we'll reuse the `orders.tsv` file from the _Terminal Basics_ tutorial (<Execute command="head orders.tsv" inline />), which contains <Link href="https://github.com/TheUpshot/chipotle/">take-out order data from Chipotle</Link>.

Time to wrangle some data!
