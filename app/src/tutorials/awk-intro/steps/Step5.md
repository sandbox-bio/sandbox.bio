<script>
import Alert from "components/Alert.svelte";
import Link from "components/Link.svelte";
import Execute from "components/Execute.svelte";
</script>

You may be surprised to learn that `awk` supports arrays (even multidimensional ones!). As an example, let's track how many of each kind of burrito was ordered in this dataset.

First, let's extract all the lines with information about burrito orders. Instead of using an `if` statement for each type of burrito (chicken, steak, carnitas, and vegetarian), we'll use the tilde operator `~` to perform a regular expression match:

<Execute command={`awk -F "\\t" '{ if($3 ~ /Burrito/) print }' orders.tsv | head`} />

Now at each line, we'll increment how often we've seen each type of burrito simply by using the syntax `arrayName[index]++`, and output the number of `Chicken Burrito`s we've encountered:

<Execute command={`awk -F "\\t" ' \\ { if($3 ~ /Burrito/) counts[$3]++ } \\ END { print counts["Chicken Burrito"] }' orders.tsv`} />

To loop through all keys and values of the array, we can simply use a `for` loop:

<Execute command={`awk -F "\\t" ' \\ { if($3 ~ /Burrito/) counts[$3]++ } \\ END { for(k in counts) print(k, counts[k]) }' orders.tsv`} />

Note that there are 6 unmarked burritos in our dataset, likely due to missing data!
