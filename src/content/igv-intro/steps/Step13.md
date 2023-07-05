<script>
import Alert from "$components/Alert.svelte";
import IGVUpdateBtn from "$components/IGVUpdateBtn.svelte";
</script>

Navigate to region `chr21:19,089,694-19,095,362`

<IGVUpdateBtn locus="chr21:19,089,694-19,095,362" />

<Alert color="primary">
	**Notes**

    * There are many reads with mismatches to the reference.
    * Read pairs are in RL pattern (instead of LR pattern).
    * Region is flanked by reads with poor mapping quality (white instead of grey).
    * There are reads with pairs on other chromosomes (coloured reads).

</Alert>
