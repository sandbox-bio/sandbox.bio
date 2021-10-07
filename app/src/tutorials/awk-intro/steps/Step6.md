<script>
import Execute from "components/Execute.svelte";
</script>

Next, let's tackle a slightly more difficult question: let's count how much money was spent on burritos, broken down by burrito filling.

Remember (or use <Execute command="head orders.tsv" inline />) that the last column contains the price. To calculate the total amount paid for each line of the file, we can multiply the `item_price` (`$NF`) by the `quantity` (`$2`) and sum those up.

But because the last field is a string, we need to remove the `$` before we can add the numbers using the `sub` command:

<Execute command={`awk -F "\\t" '{ \\ sub(/\\$/, "", $NF); \\ if($3 ~ /Burrito/) counts[$3] += $NF * $2 } \\ END { for(k in counts) print(k, counts[k]) }' orders.tsv`} />

A few notes:

* The `sub` function modifies the variable `$NF` (temporarily), and returns the number of substitutions, i.e. it does **not** return the modified string, as you might expect.
* We have to escape the `$` character since it's otherwise considered a valid regular expression.
