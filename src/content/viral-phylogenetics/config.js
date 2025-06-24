// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "viral-phylogenetics",
	name: "Viral Phylogenetics",
	icon: "virus",
	date: "July 2025",
	subtitle: `by <a href="https://www.linkedin.com/in/madison-ritter-b1374518a" target="_blank">Maddie Ritter</a> and <a href="https://www.linkedin.com/in/kyra-fetter" target="_blank">Kyra Fetter</a>`,
	description: "Perform viral phylogenetics analysis using real SARS-COV2 whole-genome sequences.",
	tags: ["MSA", "phylogeny", "FastTree", "lsd2"],
	tools: ["FastTree", "lsd2", "less", "nw_display", "ViralMSA.py"],
	difficulty: ["intermediate"],
	steps: [
		{ name: "Viral Phylogenetics", component: Intro },
		{ name: "Multiple Sequence Alignment", component: Step1, subtitle: "The sequencing reads", header: true },
		{ name: "Multiple Sequence Alignment", component: Step2, subtitle: "Run multiple sequence alignment" },
		{ name: "Generate Unrooted Tree", component: Step3, subtitle: "Use FastTree to generate a tree from alignment data", header: true },
		{ name: "Generate Rooted Tree", component: Step4, subtitle: "Use LSD2 and dates file to generate a rooted tree" },
		{ name: "Can We Date the MRCA of HIV-1?", component: Step5, subtitle: "Let's perform a viral phylogenetic analysis on HIV-1 genomes", header: true },
		{ name: "The end", component: Conclusion, subtitle: "Conclusion", header: true }
	],
	files: ["hiv1_sequences.fas", "hiv1_dates.txt", "hiv1_outgroups.txt", "hiv1_reference.fas", "sarscov2_sequences.fas", "sarscov2_dates.txt", "sarscov2_outgroup.txt", "sarscov2_reference.fas"],
};
