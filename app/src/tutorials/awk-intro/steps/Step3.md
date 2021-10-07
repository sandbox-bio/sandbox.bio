<script>
import Alert from "components/Alert.svelte";
import Link from "components/Link.svelte";
import Execute from "components/Execute.svelte";
</script>

So far, the filtering we did using `awk` could have been done using other command-line tools such as `cut` and `grep`. Let's now tackle something a bit more involved: counting the total number of times in this dataset where someone ordered a chicken bowl.

To do so, we'll use `awk`'s `BEGIN` and `END` statement, which let you run some code once before processing your file, and once after all rows in the file have been processed:

<Execute command={`awk -F "\\t" ' \\ BEGIN { sum = 0 } \\ { if($3 == "Chicken Bowl") sum += $2 } \\ END { print(sum) }' orders.tsv`} />

How this works:

* We initialize `sum` to 0
* If the customer ordered a `Chicken Bowl`, we increment `sum` by the number of bowls ordered
* At the end, print the total `sum`

Note that `awk` automatically initalizes variables for you, so `sum = 0` is not strictly necessary (but preferable for clarity):

<Execute command={`awk -F "\\t" ' \\ { if($3 == "Chicken Bowl") sum += $2 } \\ END { print(sum) }' orders.tsv`} />

As a side note, `awk` doesn't actually need a file to work on if you only provide a `BEGIN` statement with no body:

<Execute command={`awk 'BEGIN{ print(5/7) }'`} />

(that's another way to do floating-point math on the command line without using the `bc` command)
