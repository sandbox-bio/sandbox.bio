<script>
import Execute from "$components/Execute.svelte";
</script>

You may be surprised to learn that `awk` supports arrays (even multidimensional ones!). As an example, let's track how many of each kind of burrito was ordered in this dataset.

First, let's extract all the lines with information about burrito orders. Instead of using an `if` statement for each type of burrito (chicken, steak, carnitas, and vegetarian), we'll use the tilde operator `~` to perform a regular expression match:

<Execute command={`awk -F "\\t" '{ if($3 ~ /Burrito/) print }' orders.tsv | head`} />

Now we can use the array `counts` to track how many of each burrito was ordered. We'll use the syntax `counts[burritoType]` to access the counts of each burrito filling, and we can output the number of chicken burritos at the end:

<Execute command={`awk -F "\\t" ' \\ { if($3 ~ /Burrito/) counts[$3] += $2 } \\ END { print counts["Chicken Burrito"] }' orders.tsv`} />

To loop through all keys and values of the array, we can simply use a `for` loop:

<Execute command={`awk -F "\\t" ' \\ { if($3 ~ /Burrito/) counts[$3] += $2 } \\ END { for(k in counts) print(k, counts[k]) }' orders.tsv`} />

Note that there are 6 unmarked burritos in our dataset, likely due to missing data!
