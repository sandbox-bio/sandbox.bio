<script>
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

The GitHub issues API has a lot of details I don't care about. I want to select multiple fields from the returned JSON document and leave the rest behind.

The easiest way to do this is using `,` to specify multiple filters:

<Execute command={`jq ' .[].title, .[].number' issues.json`} />

But this is returning the results of one selection after the other. To change the ordering, I can factor out the array selector:

<Execute command={`jq '.[] |  .title, .number' issues.json`} />

This refactoring uses a pipe (`|`), which I'll talk about shortly, and runs my object selectors (`.title` and `.number`) on each array element.

If you wrap the query in the array constructor you get this:

<Execute command={`jq '[ .[] |  .title, .number ]' issues.json`} />

But this still isn't the JSON document I need. To get these values into a proper JSON object, I need an object constructor <code>&#123;...&#125;</code>.
