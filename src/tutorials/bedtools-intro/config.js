// Steps
import BedtoolsIntro from "./steps/Intro.md";
import BedtoolsStep1 from "./steps/Step1.md";
import BedtoolsStep2 from "./steps/Step2.md";
import BedtoolsStep3 from "./steps/Step3.md";
import BedtoolsStep4 from "./steps/Step4.md";
import BedtoolsStep5 from "./steps/Step5.md";
import BedtoolsStep6 from "./steps/Step6.md";
import BedtoolsStep7 from "./steps/Step7.md";
import BedtoolsStep8 from "./steps/Step8.md";
import BedtoolsStep9 from "./steps/Step9.md";
import BedtoolsStep10 from "./steps/Step10.md";
import BedtoolsStep11 from "./steps/Step11.md";
import BedtoolsStep12 from "./steps/Step12.md";
import BedtoolsStep13 from "./steps/Step13.md";
// Exercises
import BedtoolsExercise1 from "./exercises/Exercise1.md";
import BedtoolsExercise2 from "./exercises/Exercise2.md";
import BedtoolsExercise3 from "./exercises/Exercise3.md";
// Data files
import BedtoolsCpgBed from "./data/cpg.bed";
import BedtoolsExonsBed from "./data/exons.bed";
import BedtoolsData1Bed from "./data/fHeart-DS15839.bed";
import BedtoolsData2Bed from "./data/fHeart-DS16621.bed";
import BedtoolsData3Bed from "./data/fSkin-DS19745.bed";
import BedtoolsGwasBed from "./data/gwas.bed";
import BedtoolsChromHMMBed from "./data/hesc.chromHmm.bed";
import BedtoolsGenomeTxt from "./data/genome.txt";

export const config = {
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

	// Lessons
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
		{ name: "Exercises", component: BedtoolsExercise1, subtitle: "Find non-exons", header: true },
		{ name: "Exercises", component: BedtoolsExercise2, subtitle: "Find exons in 500kb regions" },
		{ name: "Exercises", component: BedtoolsExercise3, subtitle: "Find flanking splice sites" },
	],

	// Files needed at runtime
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
}
