<script>
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

The issues list I was looking at has many low-quality issues in it. Let's say I want to grab all the items that are labeled. This would let me skip all the drive-by fix-my-problem issues.

Unfortunately, it's impossible to do this with the GitHub API unless you specify all the possible labels in your query. However, I can easily do this query on the command line by filtering our results with `jq`. However, to do so, I'm going to need a couple more `jq` functions.

My query so far looks like this:

<Execute command={`jq '[ .[] | { title: .title, number: .number, labels: .labels | length } ]' issues.json`} />

The first thing I can do is simplify it using map.

<Execute command={`jq 'map({ title: .title, number: .number, labels: .labels | length })' issues.json`} />

`map(...)` lets you unwrap an array, apply a filter and then rewrap the results back into an array. You can think of it as a shorthand for `[ .[] | ... ]` and it comes up quite a bit in my experience, so it's worth it committing to memory.

I can combine that with a select statement that looks like this:

`map(select(.labels > 0))`

`select` is a built-in function that takes a boolean expression and only returns elements that match. It's similar to the `WHERE` clause in a SQL statement or array filter in JavaScript.

Like `map`, I find `select` comes up quite a bit, so while you may have to come back to this article or google it the first few times you need it, with luck, it will start to stick to your memory after that.

Putting this all together looks like this:

<Execute command={`jq 'map({ title: .title, number: .number, labels: .labels | length }) | \\ map(select(.labels > 0))' issues.json`} />

This uses three object indexes, two maps, two pipes, a `length` function, and a `select` predicate. But if you've followed along, this should all make sense. It's all just composing together filters until you get the result you need.

<Alert>
	**What I Learned**:

	`jq` lets you select elements by starting with a `.` and accessing keys and arrays like it's a JavaScript Object (which it is). This feature uses the Object and Array index `jq` creates of a JSON document and look like this:

	`jq '.key[0].subkey[2:3].subsubkey'`

	`jq` programs can contain object constructors `&#123; ... &#125;` and array constructors `[ ... ]`. You use these when you want to wrap back up something you've pulled out of a JSON document using the above indexes:

	`jq '[&#123; key1: .key1, key2: .key2 &#125;]'`

	`jq` contains built-in functions (`length`, `sort`, `select`, `map`) and pipes (`|`), and you can compose these together just like you can combine pipes and filters at the command line:

	`jq 'map(&#123; order-of-magitude: .items | length | tostring | length &#125;)`
</Alert>
