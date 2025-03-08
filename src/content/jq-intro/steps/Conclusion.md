<script>
import Execute from "$components/Execute.svelte";
</script>

This concludes the tutorial on using `jq`!

> **What I Learned**:
> 
> `jq` lets you select elements by starting with a `.` and accessing keys and arrays like it's a JavaScript Object (which it is). This > feature uses the Object and Array index `jq` creates of a JSON document and look like this:
> 
> ```bash
> jq '.key[0].subkey[2:3].subsubkey'
> ```
> 
> `jq` programs can contain object constructors <code>&#123; ... &#125;</code> and array constructors `[...]`. You use these when you want > to wrap back up something you've pulled out of a JSON document using the above indexes:
> 
> ```bash
> jq '[&lbrace; key1: .key1, key2: .key2 }]'
> ```
> 
> `jq` contains built-in functions (`length`, `sort`, `select`, `map`) and pipes (`|`), and you can compose these together just like you can > combine pipes and filters at the command line:
> 
> ```bash
> jq 'map(&lbrace; order-of-magitude: .items | length | tostring | length })
> ```
