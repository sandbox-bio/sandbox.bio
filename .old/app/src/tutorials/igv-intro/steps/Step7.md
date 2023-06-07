<script>
import Alert from "components/Alert.svelte";
import Link from "components/Link.svelte";
import IGVUpdateBtn from "components/IGVUpdateBtn.svelte";
</script>

#### Example 1

Navigate to position `chr21:19,518,412-19,518,497`:

<IGVUpdateBtn locus="chr21:19,518,412-19,518,497" />

Center on the A within the homopolymer run (`chr21:19,518,470`):

<IGVUpdateBtn locus="chr21:19,518,470" />

Now color by read strand and sort alignments by base at the center. Note the poor base quality support in forward reads.

#### Example 2

Center on the one base deletion (`chr21:19,518,452`):

<IGVUpdateBtn locus="chr21:19,518,452" />

Sort alignments by base, and notice the alternating insertions and deletions in reverse strand reads

<Alert color="primary">
	**Notes**

	* The alternate allele is either a deletion or an insertion of one or two Ts.
	* The remaining bases are mismatched, because the alignment is now out of sync.
	* The dpSNP entry at this location (<Link href="https://www.ncbi.nlm.nih.gov/snp/?term=rs74604068">rs74604068</Link>) is an `A->T`, and in all likelihood an artifact.
	* i.e. the common variants from dbSNP include some cases that are actually common misalignments caused by repeats.
</Alert>
