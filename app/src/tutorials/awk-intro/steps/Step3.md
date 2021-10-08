<script>
import Execute from "components/Execute.svelte";
</script>

Let's now tackle something more involved: counting the number of times someone ordered a chicken bowl.

To do so, we'll use `awk`'s `BEGIN` and `END` statements, which let you run some code once before processing your file (`BEGIN`), and once after all rows in the file have been processed (`END`):

<Execute command={`awk -F "\\t" ' \\ BEGIN { sum = 0 } \\ { if($3 == "Chicken Bowl") sum += $2 } \\ END { print(sum) }' orders.tsv`} />

Let's break this down:

* **Begin block**: We initialize the variable `sum` to 0
* **Body block**: For each line in the file, if the customer ordered a `Chicken Bowl`, we increment `sum` by the number of bowls ordered
* **End block**: Once we processed all lines, print the total `sum`
