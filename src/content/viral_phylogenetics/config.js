// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "viral-phylogenetics",
	name: "Viral Phylogenetics",
	subtitle: `by <a href="https://www.linkedin.com/in/madison-ritter-b1374518a" target="_blank">Maddie Ritter</a> and <a href="https://www.linkedin.com/in/kyra-fetter" target="_blank">Kyra Fetter</a>`,
	description: "Perfrom viral phylogenetics analysis using real HIV-1 whole-genome sequences.",
	tags: ["MSA", "phylogeny", "FastTree", "LSD2"],
	tools: ["mafft", "FastTree", "LSD2", "less", "ls", "echo", "head", "cat"],
	difficulty: ["intermediate"],
	steps: [
		{ name: "Viral Phylogenetics", component: Intro },
		{ name: "Data Analysis", component: Step1, subtitle: "Align the reads", header: true },
		{ name: "Data Analysis", component: Step2, subtitle: "Perform phylogenetic inference", header: true },
		{ name: "Data Analysis", component: Step3, subtitle: "Root and date the phylogeny" },
		{ name: "The end", component: Conclusion, subtitle: "Conclusion", header: true }
	],
	files: ["hiv1_sequences.fas", "hiv1_dates.txt"],
	//init: `REF_FASTA=reference.fasta
//PRIMER_BED=primer.bed`
};
