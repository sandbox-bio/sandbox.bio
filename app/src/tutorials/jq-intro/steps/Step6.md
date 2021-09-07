<script>
import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

The next problem I have is that I want to summarize some this JSON data. Each issue returned by GitHub has a collection of labels:

<Execute command={`jq '{ title: .title, number: .number, labels: .labels }' issue.json`} />

If I want those labels in alphabetical order I can use the built in `sort` function. It works like this:

<Execute command={`echo '["3","2","1"]' | jq 'sort'`} />

This is similar to how I would sort an array in JavaScript:

<code>
const l = ["3","2","1"];
l.sort();
</code>

Other built-ins that mirror JavaScript functionality are available, like `length`, `reverse`, and `tostring` and they can all be used in a similar way:

<Execute command={`echo '["3","2","1"]' | jq 'reverse'`} />

<Execute command={`echo '["3","2","1"]' | jq 'length'`} />

If I can combine these built-ins with the selectors I've built up so far, I'll have solved my label sorting problem. So I'll show that next.

<Alert>
	**What I Learned: Array-Index**:

	`jq` has many built-in functions. There are probably too many to remember but the built-ins tend to mirror JavaScript functions, so give those a try before heading to the <Link href="https://stedolan.github.io/jq/manual/#Builtinoperatorsandfunctions">jq manual</Link>, and you might get lucky.
</Alert>
