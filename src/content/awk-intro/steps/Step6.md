<script>
import Execute from "$components/Execute.svelte";
</script>

Next, let's tackle a slightly more difficult question: let's count how much money was spent on burritos, broken down by filling.

Remember (or use <Execute command="head orders.tsv" inline />) that the last column contains the price. Therefore, to calculate the total amount paid in each line of the file, we can multiply the `item_price` (`$NF`) by the `quantity` (`$2`) and sum those up.

But because the price column is actually a string, we can't simply sum the numbers up; we first need to remove the dollar sign (`$`) using the string replace command `sub`:

<Execute command={`awk -F "\\t" '{ \\ sub(/\\$/, "", $NF); \\ if($3 ~ /Burrito/) counts[$3] += $NF * $2 } \\ END { for(k in counts) print(k, counts[k]) }' orders.tsv`} />

A few notes:

- The `sub` function **modifies** the variable `$NF` (temporarily) and returns the number of substitutions, i.e. it does **not** return the modified string, as you might expect.
- We have to escape the `$` character since it's otherwise considered a valid regular expression.
