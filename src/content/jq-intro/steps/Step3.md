<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

Once you start using the array index to select elements, you have a new problem. The data returned won't be a valid JSON document. In the example above, the issue titles were new line delimited:

<Execute command={`jq '.[].title' issues.json`} />

In fact, whenever you ask `jq` to return an unwrapped collection of elements, it prints them each on a new line. You can see this by explicitly asking `jq` to ignore its input and instead return two numbers:

<Execute command={`echo '""' | jq '1,2'`} />

You can resolve this the same way you would turn the text `1,2` into an array in JavaScript: By wrapping it in an array constructor `[ ... ]`.

<Execute command={`echo '""' | jq '[1,2]'`} />

Similarly, to put a generated collection of results into a JSON array, you wrap it in an array constructor `[ ... ]`.

My GitHub issue title filter (`.[].title`) then becomes `[ .[].title ]` like this:

<Execute command={`jq '[ .[].title ]' issues.json`} />

Now I have a valid JSON document.

> **What I Learned: Array Constructors**:
> 
> If your `jq` query returns more than one element, they will be returned newline delimited.
> 
> <Execute command={`echo '[{"a":"b"},{"a":"c"}]' | jq '.[].a'`} />
> 
> To turn these values into a JSON array, what you do is similar to creating an array in JavaScript: You wrap the values in an array > constructor (`[...]`).
> 
> It ends up looking something like this:
> 
> <Execute command={`echo '[{"a":"b"},{"a":"c"}]' | jq '[ .[].a ]'`} />
