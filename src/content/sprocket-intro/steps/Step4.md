<script>
import Execute from "$components/Execute.svelte";
</script>

By default, Sprocket runs each task in a Docker container. This sandbox has no Docker, so `sprocket.toml` selects the **local** backend:

<Execute command="cat sprocket.toml" />

Sprocket reads `sprocket.toml` from the current directory automatically. Now run the workflow:

<Execute command="sprocket run overlaps.wdl @inputs.json" />

Sprocket runs the three `count_overlaps` calls and then `summarize`. In this emulated machine that takes about 10–15 seconds. When it finishes, it prints the workflow outputs as JSON.

> Sprocket warns that tasks won't run in a container. That's expected here. The `container` lines in the WDL are used when you run the same workflow with Docker, Apptainer, or on an HPC cluster.
