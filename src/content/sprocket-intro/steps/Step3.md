<script>
import Execute from "$components/Execute.svelte";
</script>

What inputs does the workflow need? Ask Sprocket:

<Execute command="sprocket inputs overlaps.wdl" />

We need one `features` file and an array of `samples`. Our `inputs.json` fills them in:

<Execute command="cat inputs.json" />

Let's validate these inputs against the workflow without running it:

<Execute command="sprocket validate overlaps.wdl @inputs.json" />

The `@` tells Sprocket to read inputs from a file. Validation catches wrong types and missing inputs before any task runs, which matters most when a workflow takes hours.
