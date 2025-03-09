// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Step7 from "./steps/Step7.md";
import Step8 from "./steps/Step8.md";
import Step9 from "./steps/Step9.md";
import Step10 from "./steps/Step10.md";
import Step11 from "./steps/Step11.md";
import Step12 from "./steps/Step12.md";
import Step13 from "./steps/Step13.md";
import Conclusion from "./steps/Conclusion.md";

const TRACK_HCC1143 = {
	url: "https://assets.sandbox.bio/tutorials/igv-intro/HCC1143.normal.21.19M-20M.bam",
	indexURL: "https://assets.sandbox.bio/tutorials/igv-intro/HCC1143.normal.21.19M-20M.bam.bai",
	height: 500,
	name: "HCC1143",
	viewAsPairs: true
};
const TRACK_GC = {
	url: "https://data.broadinstitute.org/igvdata/annotations/hg19/hg19.gc5base.tdf",
	name: "GC Track"
};
const TRACK_REPEATS = {
	name: "Repeat Masker (rmsk)",
	type: "annotation",
	format: "rmsk",
	displayMode: "EXPANDED",
	url: "https://s3.amazonaws.com/igv.org.genomes/hg19/rmsk.txt.gz",
	indexURL: "https://s3.amazonaws.com/igv.org.genomes/hg19/rmsk.txt.gz.tbi",
	visibilityWindow: 1000000
};

export const config = {
	id: "igv-intro",
	name: "Visualize variants with IGV",
	subtitle: `by Malachi Griffith, Sorana Morrissy, Jim Robinson, Ben Ainscough, Jason Walker, and Obi Griffith`,
	description: "Distinguish real variants from artifacts using the IGV genome browser.",
	tags: ["igv"],
	tools: [],
	difficulty: ["beginner"],
	igv: true,
	steps: [
		{ name: "Visualize variants with IGV", component: Intro },
		{ name: "Get familiar with the interface", component: Step1, subtitle: "Navigation", header: true },
		{ name: "Get familiar with the interface", component: Step2, subtitle: "Colors" },
		{ name: "Get familiar with the interface", component: Step3, subtitle: "Gene model" },
		{ name: "Get familiar with the interface", component: Step4, subtitle: "Loading Read Alignments" },
		{ name: "Get familiar with the interface", component: Step5, subtitle: "Visualizing read alignments" },
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step6, subtitle: "Two neighbouring SNPs", header: true },
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step7, subtitle: "Homopolymer region with indel" },
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step8, subtitle: "Coverage by GC" },
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step9, subtitle: "Heterozygous SNPs on different alleles" },
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step10, subtitle: "Low mapping quality" },
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step11, subtitle: "Homozygous deletion" },
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step12, subtitle: "Mis-alignment" },
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step13, subtitle: "Translocation" },
		{ name: "The End", component: Conclusion, header: true }
	],
	files: [],
	igvConfig: {
		default: {
			locus: "all",
			genome: "hg19",
			tracks: [],
			showCenterGuide: true
		},
		2: { locus: "chr1:10,000-11,000" },
		3: { locus: "BRCA1" },
		5: { tracks: [TRACK_HCC1143] },
		6: { locus: "chr21:19,480,041-19,480,386", tracks: [TRACK_HCC1143] },
		7: { locus: "chr21:19,479,301-19,479,341", tracks: [TRACK_HCC1143] },
		8: { locus: "chr21:19,518,432-19,518,472", tracks: [TRACK_HCC1143] },
		9: { locus: "chr21:19,611,925-19,631,555", tracks: [TRACK_HCC1143, TRACK_GC] },
		10: { locus: "chr21:19,666,881-19,666,921", tracks: [TRACK_HCC1143] },
		11: { locus: "chr21:19,800,320-19,818,162", tracks: [TRACK_HCC1143, TRACK_REPEATS] },
		12: { locus: "chr21:19,324,469-19,331,468", tracks: [TRACK_HCC1143, TRACK_REPEATS] },
		13: { locus: "chr21:19,102,154-19,103,108", tracks: [TRACK_HCC1143] }
	}
};
