<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

If you `curl` the GitHub Issues API, you will get back an array of issues:

<Execute command={`curl "https://api.github.com/repos/stedolan/jq/issues?per_page=5" > issues.json; cat issues.json`} />

To get a specific element in the array, give `jq` an index:

<Execute command={`jq '.[4]' issues.json`} />

<Alert>
	**Side Note: Array Indexing in `jq`**:

    Array indexing has some helpful convenience syntax.

    You can select ranges:

    <Execute command={`echo "[1,2,3,4,5]" | jq '.[2:4]'`} />

    You can select one sided ranges:

    <Execute command={`echo "[1,2,3,4,5]" | jq '.[2:]'`} />

    Also, you can use negatives to select from the end:

    <Execute command={`echo "[1,2,3,4,5]" | jq '.[-2:]'`} />

</Alert>

You can use the array index with the object index:

<Execute command={`jq '.[4].title' issues.json`} />

And you can use `[]` to get all the elements in the array. For example, here is how I would get the titles of the issues returned by my API request:

<Execute command={`jq '.[].title' issues.json`} />

<Alert>
	**What I Learned: Array-Index**:

    `jq` lets you select the whole array `[]`, a specific element `[3]`, or ranges `[2:5]` and combine these with the object index if needed.

    It ends up looking something like this:

    <code>jq '.key[].subkey[2]'</code>

</Alert>

<Alert>
	**Side Note: Removing Quotes From JQ Output**:

    The -r option in `jq` gives you raw strings if you need that.

    <Execute command={`echo '["1","2","3"]' | jq -r '.[]'`} />

    The `-j` option (for join) can combine together your output.

    <Execute command={`echo '["1","2","3"]' | jq -j '.[]'`} />

</Alert>
