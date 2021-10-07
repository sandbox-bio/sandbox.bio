<script>
import Alert from "components/Alert.svelte";
import Link from "components/Link.svelte";
import Execute from "components/Execute.svelte";
</script>

We can also use `awk` to fetch rows that match a certain pattern. For example, if we're only interested in chicken bowls, we can specify the conditions under which a row should be output:

<Execute command={`awk -F "\\t" '$3 == "Chicken Bowl"' orders.tsv | head`} />

Or to only retrieve certain columns:

<Execute command={`awk -F "\\t" '$3 == "Chicken Bowl" { print $2, $3 }' orders.tsv | head`} />

An alternative to this pattern is to use an `if` statement to decide whether to `print` a line or not (I prefer this approach since it looks more like other programming languages):

<Execute command={`awk -F "\\t" '{ if($3 == "Chicken Bowl") print $2, $3 }' orders.tsv | head`} />
