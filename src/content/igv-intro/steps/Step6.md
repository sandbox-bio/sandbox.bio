<script>
import Alert from "$components/Alert.svelte";
import IGVUpdateBtn from "$components/IGVUpdateBtn.svelte";
</script>

In this section we will be looking in detail at 8 positions in the genome, and determining whether they represent real events or artifacts.

##### Two neighbouring SNPs

Navigate to region `chr21:19,479,237-19,479,814`:

<IGVUpdateBtn locus="21:19,479,237-19,479,814" />

Note two heterozygous variants, one corresponds to a known dbSNP (`G/T` on the right) the other does not (`C/T` on the left).

Zoom in and center on the `C/T` SNV on the left:

<IGVUpdateBtn locus="21:19,479,321" />

Sort alignments by `base`: right-click on a center base and choose `Sort by base`.

Color alignments by `read strand`: click the track gear icon on the right and choose `Color by read strand`.

<Alert color="primary">
	**Notes**

	* Reads have high base quality scores except one (where the alt allele is the last base of the read). Note that `igv.js` dims the color of a base if it has low quality.
	* Reads have good mapping qualities, show no strand bias, and the allele frequency is consistent with a heterozygous mutation.
</Alert>

<Alert color="info">
	**Question**

	How does `Color by read strand` help?
</Alert>
