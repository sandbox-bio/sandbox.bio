<script>
import Alert from "$components/Alert.svelte";
import IGVUpdateBtn from "$components/igv/IGVUpdateBtn.svelte";
</script>

Navigate to region `chr21:19,666,881-19,666,921`:

<IGVUpdateBtn locus="chr21:19,666,881-19,666,921" />

Sort by base at the center position.

There are two SNPs in this region, `A/T` on the left and `AG` on the right.

<Alert color="primary">
	**Note**

    There is no linkage between alleles for these two SNPs because reads covering both only contain one or the other

</Alert>
