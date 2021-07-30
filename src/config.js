// TODO: move this so that each tutorial has its own file!
// Introduction to bedtools
import BedtoolsIntro from "./tutorials/bedtools-intro/Intro.svelte";
import BedtoolsStep1 from "./tutorials/bedtools-intro/Step1.svelte";
import BedtoolsStep2 from "./tutorials/bedtools-intro/Step2.svelte";
import BedtoolsStep3 from "./tutorials/bedtools-intro/Step3.svelte";
import BedtoolsStep4 from "./tutorials/bedtools-intro/Step4.svelte";
import BedtoolsStep5 from "./tutorials/bedtools-intro/Step5.svelte";
import BedtoolsStep6 from "./tutorials/bedtools-intro/Step6.svelte";
import BedtoolsStep7 from "./tutorials/bedtools-intro/Step7.svelte";
import BedtoolsStep8 from "./tutorials/bedtools-intro/Step8.svelte";
import BedtoolsStep9 from "./tutorials/bedtools-intro/Step9.svelte";
import BedtoolsStep10 from "./tutorials/bedtools-intro/Step10.svelte";
import BedtoolsStep11 from "./tutorials/bedtools-intro/Step11.svelte";
import BedtoolsStep12 from "./tutorials/bedtools-intro/Step12.svelte";
import BedtoolsStep13 from "./tutorials/bedtools-intro/Step13.svelte";
import BedtoolsConclusion from "./tutorials/bedtools-intro/Conclusion.svelte";
// Bedtools data files
import BedtoolsCpgBed from "./tutorials/bedtools-intro/data/cpg.bed";
import BedtoolsExonsBed from "./tutorials/bedtools-intro/data/exons.bed";
import BedtoolsData1Bed from "./tutorials/bedtools-intro/data/fHeart-DS15839.bed";
import BedtoolsData2Bed from "./tutorials/bedtools-intro/data/fHeart-DS16621.bed";
import BedtoolsData3Bed from "./tutorials/bedtools-intro/data/fSkin-DS19745.bed";
import BedtoolsGwasBed from "./tutorials/bedtools-intro/data/gwas.bed";
import BedtoolsChromHMMBed from "./tutorials/bedtools-intro/data/hesc.chromHmm.bed";
import BedtoolsGenomeTxt from "./tutorials/bedtools-intro/data/genome.txt";

export const config = {
	"tutorials": [
		{
			// Metadata
			id: "bedtools-intro",
			name: "Introduction to bedtools",
			description: "Explore, analyze, and manipulate genomic interval <code>.bed</code> files.",
			tools: ["bedtools"],
			difficulty: ["beginner"],
			adapted_from: {
				"name": "quinlanlab.org/tutorials/bedtools",
				"link": "http://quinlanlab.org/tutorials/bedtools/bedtools.html"
			},
			// Tutorial
			steps: [
				{ name: "Introduction to bedtools", component: BedtoolsIntro },
				{ name: "The data", component: BedtoolsStep1 },
				{ name: "Bedtools help", component: BedtoolsStep2 },
				{ name: "bedtools intersect", component: BedtoolsStep3, subtitle: "Intersect regions", header: true },
				{ name: "bedtools intersect", component: BedtoolsStep4, subtitle: "Find overlaps" },
				{ name: "bedtools intersect", component: BedtoolsStep5, subtitle: "Count overlaps" },
				{ name: "bedtools intersect", component: BedtoolsStep6, subtitle: "Faster with sorted data" },
				{ name: "bedtools intersect", component: BedtoolsStep7, subtitle: "Intersect multiple files" },
				{ name: "bedtools merge", component: BedtoolsStep8, subtitle: "Merge regions", header: true },
				{ name: "bedtools merge", component: BedtoolsStep9, subtitle: "Merge and count intervals" },
				{ name: "bedtools merge", component: BedtoolsStep10, subtitle: "Merge nearby features" },
				{ name: "bedtools complement", component: BedtoolsStep11, subtitle: "Find uncovered regions", header: true},
				{ name: "bedtools genomecov", component: BedtoolsStep12, subtitle: "Measure genome-wide coverage", header: true },
				{ name: "bedtools jaccard", component: BedtoolsStep13, subtitle: "Measure dataset similarity", header: true },
				{ name: "The end", component: BedtoolsConclusion, header: true },
			],
			// Files needed
			files: [
				{ name: "cpg.bed", contents: BedtoolsCpgBed },
				{ name: "exons.bed", contents: BedtoolsExonsBed },
				{ name: "fHeart-DS15839.bed", contents: BedtoolsData1Bed },
				{ name: "fHeart-DS16621.bed", contents: BedtoolsData2Bed },
				{ name: "fSkin-DS19745.bed", contents: BedtoolsData3Bed },
				{ name: "gwas.bed", contents: BedtoolsGwasBed },
				{ name: "hesc.chromHmm.bed", contents: BedtoolsChromHMMBed },
				{ name: "genome.txt", contents: BedtoolsGenomeTxt }
			]
		},
		{
			"id": "bowtie2-intro",
			"name": "Introduction to bowtie2",
			"description": "Align DNA sequencing reads from <code>.fastq</code> files to the Lambda phage reference genome."
		}
	]
}
