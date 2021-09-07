<script>
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

I'm going to use `jq` to filter the data returned by the GitHub repository API:

<Execute command={`curl "https://api.github.com/repos/stedolan/jq" > repo.json`} />

`jq` lets us treat the JSON document as an object and select elements inside of it.

Here is how I filter the JSON document to select the value of the `name` key:

<Execute command={`jq '.name' repo.json`} />

Similarly, for selecting the value of the `owner` key:

<Execute command={`jq '.owner' repo.json`} />

You can drill in as far as you want like this:

<Execute command={`jq '.owner.login' repo.json`} />

<Alert color="secondary">
	**What I Learned: Object Identifier-Index**:
	
	`jq` lets you select elements in a JSON document like it's a JavaScript object. Just start with `.` (for the whole document) and drill down to the value you want. It ends up looking something like this:
	
	<code>jq '.key.subkey.subsubkey'</code>
</Alert>
