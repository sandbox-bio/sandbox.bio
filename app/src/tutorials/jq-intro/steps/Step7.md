<script>
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

Before I can use `sort` to sort the labels from my GitHub API request, I need to explain how pipes and filters work in `jq`.

`jq` is a filter in the UNIX command line sense. You pipe (`|`) a JSON document to it, and it filters it and outputs it to standard out. I could easily use this feature to chain together `jq` invocations like this:

<Execute command={`echo '{"title": "JQ Select"}' | \\ jq '.title | length'`} />

Here are some more examples:

* `.title | length` will return the length of the title
* `.number | tostring` will return the issue number as a string
* `.[] | .key` will return the values of key key in the array (this is equivalent to this `.[].key`)

This means that sorting my labels array is simple. I can just change `.labels` to `.labels | sort`:

<Execute command={`jq '{ title: .title, number: .number, labels: .labels | sort }' issue.json`} />

And if you want just a label count that is easy as well:

<Execute command={`jq '{ title: .title, number: .number, labels: .labels | length }' issue.json`} />

<Alert>
	**What I Learned: Pipes and Filters**:

	Everything in `jq` is a filter that you can combine with pipes (`|`). This mimics the behavior of a UNIX shell.

	You can use the pipes and the `jq` built-ins to build complicated transformations from simple operations.

	It ends up looking something like this:

	* `jq '.key1.subkey2[] | sort'`
	* `jq '.key2.subkey | length'`
	* `jq '.key3 | floor | tostring | length'`
</Alert>
