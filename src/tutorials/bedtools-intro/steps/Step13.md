<script>
import Execute from "../../Execute.svelte";
</script>


	We will use the bedtools implementation of a Jaccard statistic to meaure the similarity of two 
	datasets. Briefly, the Jaccard statistic measures the ratio of the number of *intersecting* base 
	pairs to the *total* number of base pairs in the two sets.  As such, the score ranges from 0.0 to 1.
	0; lower values reflect lower similarity, whereas higher values reflect higher similarity.

	Let's walk through an example: we would expect the Dnase hypersensivity sites to be rather similar 
	between two samples of the **same** fetal tissue type.  Let's test:


<Execute command={"bedtools jaccard \\ -a fHeart-DS16621.bed \\ -b fHeart-DS15839.bed"} />




	But what about the similarity of two <strong>different</strong> tissue types?


<Execute command={"bedtools jaccard \\ -a fHeart-DS16621.bed \\ -b fSkin-DS19745.bed"} />




	Hopefully this demonstrates how the Jaccard statistic can be used as a simple statistic to reduce the dimensionality of the comparison between two large (e.g., often containing thousands or millions of intervals) feature sets.

