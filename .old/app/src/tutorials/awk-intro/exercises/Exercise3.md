<script>
// Solution:
//    awk -F "\t" '$3 ~ /Burrito/{ if($4 ~ /Guacamole/) with[$3] += $2; else without[$3] += $2; } END { for(k in with) print(k, with[k] / (with[k] + without[k])) }' orders.tsv > burritos_guac.tsv

import Exercise from "components/Exercise.svelte";

let criteria = [
{
	name: "File <code>burritos_guac.tsv</code> exists",
	checks: [{
		type: "file",
		path: "burritos_guac.tsv",
		action: "exists"
	}]
},
{
	name: "File <code>burritos_guac.tsv</code> contains, on each line, the burrito type, followed by the percentage of those burritos whose order include guacamole (0 to 1)",
	checks: [{
		type: "file",
		path: "burritos_guac.tsv",
		action: "contents",
		commandExpected: `awk -F "\t" '$3 ~ /Burrito/{ if($4 ~ /Guacamole/) with[$3] += $2; else without[$3] += $2; } END { for(k in with) print(k, with[k] / (with[k] + without[k])) }' orders.tsv`
	}]
}];
</script>

For each type of burrito, output the percentage (between 0 and 1) of burritos ordered that included guacamole. In those cases when a burrito contained guacamole, column `$4` will contain the string `Guacamole`.

<Exercise {criteria} />
