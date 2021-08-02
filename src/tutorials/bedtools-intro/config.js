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
// Exercises
import Exercise1 from "./exercises/Exercise1.md";
import Exercise2 from "./exercises/Exercise2.md";
import Exercise3 from "./exercises/Exercise3.md";
import Exercise4 from "./exercises/Exercise4.md";
// Data files
import CpgBed from "./data/cpg.bed";
import ExonsBed from "./data/exons.bed";
import Data1Bed from "./data/fHeart-DS15839.bed";
import Data2Bed from "./data/fHeart-DS16621.bed";
import Data3Bed from "./data/fSkin-DS19745.bed";
import GwasBed from "./data/gwas.bed";
import ChromHMMBed from "./data/hesc.chromHmm.bed";
import GenomeTxt from "./data/genome.txt";

export const config = {
	// Metadata
	id: "bedtools-intro",
	name: "Introduction to bedtools",
	description: "Explore, analyze, and manipulate genomic interval <code>.bed</code> files.",
	tools: ["bedtools"],
	difficulty: ["beginner"],

	// Lessons
	steps: [
		{ name: "Introduction to bedtools", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "Bedtools help", component: Step2 },
		{ name: "bedtools intersect", component: Step3, subtitle: "Intersect regions", header: true },
		{ name: "bedtools intersect", component: Step4, subtitle: "Find overlaps" },
		{ name: "bedtools intersect", component: Step5, subtitle: "Count overlaps" },
		{ name: "bedtools intersect", component: Step6, subtitle: "Faster with sorted data" },
		{ name: "bedtools intersect", component: Step7, subtitle: "Intersect multiple files" },
		{ name: "bedtools merge", component: Step8, subtitle: "Merge regions", header: true },
		{ name: "bedtools merge", component: Step9, subtitle: "Merge and count intervals" },
		{ name: "bedtools merge", component: Step10, subtitle: "Merge nearby features" },
		{ name: "bedtools complement", component: Step11, subtitle: "Find uncovered regions", header: true},
		{ name: "bedtools genomecov", component: Step12, subtitle: "Measure genome-wide coverage", header: true },
		{ name: "bedtools jaccard", component: Step13, subtitle: "Measure dataset similarity", header: true },
		{ name: "Exercises", component: Exercise1, subtitle: "Find non-exons", header: true },
		{ name: "Exercises", component: Exercise2, subtitle: "Find exons in 500kb regions" },
		{ name: "Exercises", component: Exercise3, subtitle: "Find flanking splice sites" },
		{ name: "Exercises", component: Exercise4, subtitle: "Calculate Jaccard statistics" }
	],

	// Files needed at runtime
	files: [
		{ name: "cpg.bed", contents: CpgBed },
		{ name: "exons.bed", contents: ExonsBed },
		{ name: "fHeart-DS15839.bed", contents: Data1Bed },
		{ name: "fHeart-DS16621.bed", contents: Data2Bed },
		{ name: "fSkin-DS19745.bed", contents: Data3Bed },
		{ name: "gwas.bed", contents: GwasBed },
		{ name: "hesc.chromHmm.bed", contents: ChromHMMBed },
		{ name: "genome.txt", contents: GenomeTxt }
	]
}
