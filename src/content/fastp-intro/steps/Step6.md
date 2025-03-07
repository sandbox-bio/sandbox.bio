<script>
import Execute from "$components/Execute.svelte";
</script>

One common FASTQ filtering step is to remove reads that are too short to be useful in further analysis (e.g. reads that are too short might not have enough information to be mapped accurately to a reference genome).

In the previous step, notice this line in the output of `fastp`:

```js
reads failed due to too short: 0
```

This means that _after_ filtering, there were no reads that were deemed too short by `fastp`. By default, the threshold for a read being too short is set to 15bp, but this can be adjusted.

For example, to only keep reads that are longer than 50bp, let's use the `--length_required` parameter:

<Execute command="fastp \ --in1 HG004_R1.fastq.gz \ --in2 HG004_R2.fastq.gz \ --length_required 50" />

With a stricter threshold, we now see that there are more reads that were discarded:

```js
reads failed due to too short: 372
```
