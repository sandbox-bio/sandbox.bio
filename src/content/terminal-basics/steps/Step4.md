<script>
import Execute from "$components/Execute.svelte";
</script>

Next, let's cover **piping**, which helps you string together many command-line tools, such that the output of one is the input of the other.

For example, to find all orders that aren't burritos, and only display the last 3, we can "pipe" the output of `grep` to `tail` using the pipe (`|`) operator:

<Execute command='grep -v "Burrito" orders.tsv | tail -n 3' />

We can also use the `wc -l` command to only **count** the number of lines without displaying them. For example, we can compare the number of Chicken vs. Steak burrito orders:

<Execute command='grep "Chicken Burrito" orders.tsv | wc -l' />

<Execute command='grep "Steak Burrito" orders.tsv | wc -l' />
