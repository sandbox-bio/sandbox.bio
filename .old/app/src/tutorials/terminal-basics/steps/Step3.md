<script>
import Execute from "components/Execute.svelte";
</script>

To generate a list of all burrito orders, we can use `grep` to filter lines using a pattern of interest:

<Execute command='grep "Burrito" orders.tsv' />

Note that `grep` is case-sensitive, i.e. the following won't return any results:

<Execute command='grep "burrito" orders.tsv' />

But we can ask `grep` to ignore the case of our pattern:

<Execute command='grep -i "burrito" orders.tsv' />

Another common pattern of inquiry is to find lines that do **not** match a pattern. To find all orders that are **not** burritos:

<Execute command='grep -v "Burrito" orders.tsv' />
