<script>
import Execute from "$components/Execute.svelte";
</script>

Each run gets its own directory under `out/runs/<workflow>/`, and `_latest` always points to the most recent one:

<Execute command="ls out/runs/peak_overlaps/_latest" />

`outputs.json` holds the workflow outputs, `inputs.json` records exactly what the run used, and `calls/` has one directory per task call:

<Execute command="ls out/runs/peak_overlaps/_latest/calls" />

Let's look at the summary table:

<Execute command="column -t out/runs/peak_overlaps/_latest/calls/summarize/attempts/0/work/summary.tsv" />

Each call keeps its own working directory, command, and logs. For example, here are the exon-overlapping peaks for the first sample:

<Execute command="head -5 out/runs/peak_overlaps/_latest/calls/count_overlaps-0/attempts/0/work/overlapping.bed" />

And the exact shell script Sprocket ran for that call:

<Execute command="cat out/runs/peak_overlaps/_latest/calls/count_overlaps-0/attempts/0/command" />
