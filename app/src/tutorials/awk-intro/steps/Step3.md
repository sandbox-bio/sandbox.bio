<script>
import Alert from "components/Alert.svelte";
import Link from "components/Link.svelte";
import Execute from "components/Execute.svelte";
</script>

So far, the filtering we did using `awk` could have been done using other command-line tools such as `cut` and `grep`. Let's now tackle something a bit more involved: counting the total number of times in this dataset where someone ordered a chicken bowl.

To do so, we'll use `awk`'s `BEGIN` and `END` statement, which let you run some code once before processing your file, and once after all rows in the file have been processed:

<Execute command={`awk -F "\\t" ' \\ BEGIN { sum = 0 } \\ { if($3 == "Chicken Bowl") sum += $2 } \\ END { print(sum) }' orders.tsv | head`} />
