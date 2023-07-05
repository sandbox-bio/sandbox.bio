<script>
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import IGVUpdateBtn from "$components/IGVUpdateBtn.svelte";
</script>

We will be using publicly available Illumina sequence data from the HCC1143 cell line. For speed, only a small portion of `chr21` will be loaded (19M:20M).

<IGVUpdateBtn
locus="21:19000000-20000000"
loadTrack={{
		url: "https://assets.sandbox.bio/tutorials/igv-intro/HCC1143.normal.21.19M-20M.bam",
		indexURL: "https://assets.sandbox.bio/tutorials/igv-intro/HCC1143.normal.21.19M-20M.bam.bai",
		name: "HCC1143"
	}}>
Load HCC1143 data into IGV
</IGVUpdateBtn>

<Alert color="primary">
	In IGV Desktop, you can load files by choosing `File` > `Load from file`.
</Alert>

<Alert>
	The HCC1143 cell line was generated from a 52 year old caucasian woman with breast cancer.
	
	Additional information on this cell line can be found here: <Link href="https://www.atcc.org/products/all/CRL-2321.aspx">HCC1143</Link> (tumor, TNM stage IIA, grade 3, primary ductal carcinoma) and <Link href="https://www.atcc.org/products/all/CRL-2362.aspx">HCC1143/BL</Link> (matched normal EBV transformed lymphoblast cell line).
</Alert>
