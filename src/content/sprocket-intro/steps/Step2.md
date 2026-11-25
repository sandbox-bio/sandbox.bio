<script>
import Execute from "$components/Execute.svelte";
</script>

Before running anything, Sprocket can catch mistakes. `sprocket check` parses the document and checks types, names and call inputs:

<Execute command="sprocket check overlaps.wdl" />

No output means no errors. To see what an error looks like, rename a variable with a typo and check again:

<Execute command="sed 's/peaks = sample/peaks = sampel/' overlaps.wdl > broken.wdl && sprocket check broken.wdl" />

Sprocket points to the exact line and name that's wrong.

`sprocket lint` also checks style and best practices:

<Execute command="sprocket lint overlaps.wdl" />

The notes say our container tags can change over time. To learn more about a rule, use `sprocket explain`:

<Execute command="sprocket explain ContainerUri" />

Finally, `sprocket format` checks that the file follows the standard WDL style:

<Execute command="sprocket format check overlaps.wdl" />
