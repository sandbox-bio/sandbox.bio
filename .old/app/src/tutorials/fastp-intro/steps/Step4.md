<script>
import Alert from "components/Alert.svelte";
import Link from "components/Link.svelte";
import Execute from "components/Execute.svelte";
</script>

`fastp` also outputs a `fastp.json` report:

<Execute command="head -n20 fastp.json" />

Let's say you're running `fastp` as part of an analysis pipeline, and you want to make an automated decision about what kind of analysis to run depending on data quality. In that scenario, it's a lot easier to parse this JSON file programmatically than any of the other outputs we've seen so far.

For example, to extract summary statistics from the report, we can use the command-line tool `jq`:

<Execute command="jq '.summary' fastp.json" />

Or to extract the number of reads that were filtered out because they had too many `N`'s in their sequences:

<Execute command="jq '.filtering_result.too_many_N_reads' fastp.json" />

<Alert>Check out <Link href="/tutorials?id=jq-intro">our `jq` tutorial</Link> for an in-depth introduction to parsing JSON files.</Alert>
