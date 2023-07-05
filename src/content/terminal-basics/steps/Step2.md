<script>
import Execute from "$components/Execute.svelte";
</script>

To view the first 10 lines of our orders dataset, you can use the `head` command:

<Execute command="head orders.tsv" />

To view a custom number of lines, use the `-n` parameter:

<Execute command="head -n 3 orders.tsv" />

From the output of `head`, you should see that the columns of this file contain the order ID, the quantity of each item, and the item name.

If we're interested in the **last** `N` lines of a file, we can instead use the aptly-named `tail` command:

<Execute command="tail orders.tsv" />

We can see that our dataset contains information about `1834` restaurant orders.

Next, we'll calculate statistics about the frequency of burrito orders!
