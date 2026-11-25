<script>
import Execute from "$components/Execute.svelte";
</script>

You don't need an input file. You can also pass inputs as `key=value` pairs on the command line. Extra values after an array input are added to that array.

Let's ask a different question: how many peaks in each sample overlap **CpG islands** instead of exons?

<Execute command="sprocket run overlaps.wdl peak_overlaps.features=cpg.bed peak_overlaps.samples=fHeart-DS15839.bed fHeart-DS16621.bed fSkin-DS19745.bed" />

<Execute command="column -t out/runs/peak_overlaps/_latest/calls/summarize/attempts/0/work/summary.tsv" />

Both runs are still on disk, each with its own `inputs.json`:

<Execute command="ls out/runs/peak_overlaps" />
