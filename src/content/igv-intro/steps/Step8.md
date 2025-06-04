<script>
import IGVUpdateBtn from "$components/igv/IGVUpdateBtn.svelte";
</script>

Navigate to position `chr21:19,611,925-19,631,555`:

<IGVUpdateBtn locus="chr21:19,611,925-19,631,555" />

Note that the range contains areas where coverage drops to zero in a few places.

Enable the collapsed view: Click the track's gear icon on the right and choose the `Display mode: squish`.

Next, we'll load a track that contains information about GC content across the genome (in IGV Desktop, this would be under `File` > `Load from server`)

<IGVUpdateBtn loadTrack={{
	url: "https://data.broadinstitute.org/igvdata/annotations/hg19/hg19.gc5base.tdf",
	name: "GC Track"
}}>
Load GC content track
</IGVUpdateBtn>

Note the concordance of coverage with GC content.

> **Question**
>
> Why are there blue and red reads throughout the alignments?
