<script>
import Execute from "$components/Execute.svelte";
</script>

##### A Little Example

Let's look at some simple examples before showing how my GitHub query can use an object constructor.

I have an array that contains my name (`["Adam", "Gordon", "Bell"]`), and I want to turn it into a JSON object like this:

```json
{
  "first": "Adam",
  "last": "Bell"
}
```

I can select the elements I need using array indexing like this:

<Execute command={`echo '["Adam", "Gordon", "Bell"]' | \\ jq -r '.[0], .[2]'`} />

To wrap those values into the shape I need, I can replace the values with the array indexes that return them:

<Execute command={`echo '["Adam", "Gordon", "Bell"]' | \\ jq '{"first":.[0], "last":.[2]}'`} />

This syntax is the same syntax for creating an object in a JSON document. The only difference is you can use the object and array queries you've built up as the values.

##### Back to GitHub

Returning to my GitHub API problem, to wrap the number and the title up into an array I use the object constructor like this:

<Execute command={`jq '[ .[] | { title: .title, number: .number}]' issues.json`} />

> **What I Learned: Object Constructors**:
>
> To put the elements you've selected back into a JSON document, you can wrap them in an object constructor.
> 
> If you were building up a JSON object out of several selectors, it would end up looking something like this:
> 
> ```shell
> jq '&lbrace; "key1": [jq filter], "key2": [jq filter] }'
> ```
> 
> Which is the same syntax for an object in a JSON document, except with jq you can use filters as values.
