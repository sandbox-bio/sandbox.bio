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
// import Conclusion from "./steps/Conclusion.md";

const TRACKS_HCC1143 = [{
	url: "https://rnabio.org/assets/module_2/HCC1143.normal.21.19M-20M.bam",
	indexURL: "https://rnabio.org/assets/module_2/HCC1143.normal.21.19M-20M.bam.bai",
	maxHeight: 500
}];

export const config = {
	id: "igv-intro",
	pwd: "igv-intro",
	name: "Visualize variants with IGV",
	subtitle: `by <a href="https://rnabio.org" target="_blank">The Griffith Lab</a>`,
	description: "Distinguish variants from artifacts.",
	tags: ["igv"],
	tools: [],
	difficulty: ["beginner"],
	igv: true,
	igvConfig: {
		default: {
			locus: "all",
			genome: "hg19",
			tracks: [],
			showCenterGuide: true,
		},
		3: {
			locus: "BRCA1"
		},
		2: {
			locus: "chr1:10,000-11,000"
		},
		5: {
			tracks: TRACKS_HCC1143
		},
		6: {
			tracks: TRACKS_HCC1143,
			locus: "21:19,480,041-19,480,386",
		},
		7: {
			tracks: TRACKS_HCC1143,
			locus: "chr21:19,479,301-19,479,341",
		},
		8: {
			tracks: TRACKS_HCC1143,
			locus: "chr21:19,518,432-19,518,472",
		},
		8: {
			tracks: TRACKS_HCC1143,
			locus: "chr21:19,518,432-19,518,472",
		},
	},
	steps: [
		{ name: "Visualize variants with IGV", component: Intro },
		{ name: "Get familiar with the interface", component: Step1, subtitle: "Navigation", header: true },
		{ name: "Get familiar with the interface", component: Step2, subtitle: "Colors" },
		{ name: "Get familiar with the interface", component: Step3, subtitle: "Gene model" },
		{ name: "Get familiar with the interface", component: Step4, subtitle: "Loading Read Alignments" },
		{ name: "Get familiar with the interface", component: Step5, subtitle: "Visualizing read alignments" },
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step6, subtitle: "Two neighbouring SNPs", header: true },
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step7, subtitle: "Homopolymer region with indel"},
		{ name: "Inspecting SNPs, SNVs, and SVs", component: Step8, subtitle: "Coverage by GC"},
	],
	files: [],
};
