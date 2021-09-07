<script>
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

##### A Little Example

Let's look at some simple examples before showing how my GitHub query can use an object constructor.

I have an array that contains my name (`["Adam", "Gordon", "Bell"]`), and I want to turn it into a JSON object like this:

<code>
&#123;
  "first": "Adam",
  "last": "Bell"
&#125;
</code>

I can select the elements I need using array indexing like this:

<Execute command={`echo '["Adam", "Gordon", "Bell"]' | \\ jq -r '.[0], .[2]'`} />

To wrap those values into the shape I need, I can replace the values with the array indexes that return them:

<Execute command={`echo '["Adam", "Gordon", "Bell"]' | \\ jq '{"first":.[0], "last":.[2]}'`} />

This syntax is the same syntax for creating an object in a JSON document. The only difference is you can use the object and array queries you've built up as the values.

##### Back to GitHub

Returning to my GitHub API problem, to wrap the number and the title up into an array I use the object constructor like this:

<Execute command={`jq '[ .[] | { title: .title, number: .number}]' issues.json`} />

<Alert>
	**What I Learned: Object Constructors**:

	To put the elements you've selected back into a JSON document, you can wrap them in an object constructor <code>&#123; ... &#125;</code>.

	If you were building up a JSON object out of several selectors, it would end up looking something like this:

	<code>jq '&#123; "key1": &lt;&lt;jq filter&gt;&gt;, "key2": &lt;&lt;jq filter&gt;&gt; &#125;'</code><br /><br />

	Which is the same syntax for an object in a JSON document, except with jq you can use filters as values.
</Alert>
