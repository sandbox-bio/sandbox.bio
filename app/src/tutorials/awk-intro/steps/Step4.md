<script>
import Execute from "components/Execute.svelte";
</script>

Note that `awk` automatically initalizes variables for you, so `sum = 0` is not strictly necessary (but preferable for clarity):

<Execute command={`awk -F "\\t" ' \\ { if($3 == "Chicken Bowl") sum += $2 } \\ END { print(sum) }' orders.tsv`} />

As a side note, `awk` doesn't actually need a file to work on if you only provide it with a `BEGIN` statement and no body:

<Execute command={`awk 'BEGIN{ print(5/7) }'`} />

(that's another way to do floating-point math on the command line without using the `bc` command!)
