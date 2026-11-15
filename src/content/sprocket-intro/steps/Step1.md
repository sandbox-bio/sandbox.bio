<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

Let's see what's in our directory: <Execute command="ls" inline />

- `overlaps.wdl` is the workflow
- `inputs.json` lists the inputs we'll run it with
- `sprocket.toml` configures Sprocket
- The `.bed` files are genomic intervals: exons, CpG islands, and DNase I peaks from three fetal tissue samples (`fHeart-*`, `fSkin-*`)

> The `.bed` files are the same subsampled files used in the <Link href="/tutorials/bedtools-intro">bedtools tutorial</Link>, which come from the <Link href="http://quinlanlab.org/tutorials/bedtools.html">Quinlan Lab bedtools tutorial</Link>. The DNase I peaks are from <Link href="https://science.sciencemag.org/content/337/6099/1190">Maurano et al. (2012)</Link>.

First, confirm Sprocket is installed:

<Execute command="sprocket --version" />

Now let's read the workflow:

<Execute command="cat overlaps.wdl" />

The workflow, `peak_overlaps`, has two tasks:

- `count_overlaps` runs `bedtools intersect` on **one** sample and writes a one-line count
- `summarize` combines those counts into a table

The `scatter` block calls `count_overlaps` once for each sample in `samples`. The calls are independent, so Sprocket can run them in parallel. Then `summarize` gathers all of their outputs.
