<script>
import Alert from "$components/Alert.svelte";
import IGVUpdateBtn from "$components/IGVUpdateBtn.svelte";
</script>

Navigate to region `chr21:19,800,320-19,818,162`:

<IGVUpdateBtn locus="chr21:19,800,320-19,818,162" />

Load the track with information about repeats across the genome (in IGV Desktop, this would be under `File` > `Load from server`):

<IGVUpdateBtn loadTrack={{
	"name": "Repeat Masker (rmsk)",
	"type": "annotation",
	"format": "rmsk",
	"displayMode": "EXPANDED",
	"url": "https://s3.amazonaws.com/igv.org.genomes/hg19/rmsk.txt.gz",
	"indexURL": "https://s3.amazonaws.com/igv.org.genomes/hg19/rmsk.txt.gz.tbi",
	"visibilityWindow": 1000000
}}>
	Load Repeat Masker track
</IGVUpdateBtn>

<Alert color="primary">
	**Note**

	Mapping quality plunges in all reads (white instead of grey). Once we load repeat elements, we see that there are two LINE elements that cause this.
</Alert>
