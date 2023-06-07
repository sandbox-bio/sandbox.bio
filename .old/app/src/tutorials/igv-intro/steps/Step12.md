<script>
import Alert from "components/Alert.svelte";
import IGVUpdateBtn from "components/IGVUpdateBtn.svelte";
</script>

Navigate to region `chr21:19,102,154-19,103,108`:

<IGVUpdateBtn locus="chr21:19,102,154-19,103,108" />

<Alert color="primary">
	**Notes**

	* This is a position where an AluY element causes mis-alignment.
	* Misaligned reads have mismatches to the reference, and well-aligned reads have partners on other chromosomes where additional ALuY elements are encoded.
	* Zoom out until you can clearly see the contrast between the difficult alignment region (corresponding to an AluY) and regions with clean alignments on either side.
</Alert>
