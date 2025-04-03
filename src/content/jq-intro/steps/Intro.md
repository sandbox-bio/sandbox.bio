<script>
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

<Alert>
	This tutorial is an interactive version of the <Link href="https://earthly.dev/blog/jq-select/">jq tutorial</Link> developed by <Link href="https://adamgordonbell.com/">Adam Gordon Bell</Link>.
</Alert>

`jq` is a lightweight, command-line JSON processor. To use it, you construct one or more filters, and it applies those filters to a JSON document.

The simplest filter is the **identity filter** which returns all its input (`.`):

<Execute command={`echo '{"key1": {"key2":"value1"}}' | jq '.'`} />

This filter is handy for just pretty-printing a JSON document. I'm going to ignore the pretty-printing and jump right into using `jq` to transform JSON documents.
