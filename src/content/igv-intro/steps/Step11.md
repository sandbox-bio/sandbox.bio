<script>
import Alert from "$components/Alert.svelte";
import IGVUpdateBtn from "$components/IGVUpdateBtn.svelte";
</script>

Navigate to region `chr21:19,324,469-19,331,468`:

<IGVUpdateBtn locus="chr21:19,324,469-19,331,468" />

Update the track settings by clicking the gear icon on the right:

* Turn on `View as pairs`
* Select `Display mode: Expand`
* Color by `Pair orientation and insert size`

Click on a red read pair to pull up information about alignments.

<Alert color="primary">
	**Notes**

	* The typical insert size of a grey read pair in the vicinity is 350bp
	* The insert size of red read pairs is 2,875bp
	* This corresponds to a homozygous deletion of 2.5kb
</Alert>
