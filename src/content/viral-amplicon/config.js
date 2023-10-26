// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "viral-amplicon",
	name: "Viral Amplicon Sequencing",
	subtitle: `by <a href="https://niema.net/" target="_blank">Niema Moshiri</a>`,
	description: "Analyze viral amplicon sequence data using a real SARS-CoV-2 dataset.",
	tags: ["ViralConsensus", "minimap", "samtools"],
	tools: ["viral_consensus", "minimap2", "samtools", "ls", "echo", "head", "cat"],
	difficulty: ["intermediate"],
	steps: [
		{ name: "Viral Amplicon Sequencing", component: Intro },
		{ name: "The Data", component: Step1 },
		{ name: "Data Analysis", component: Step2, subtitle: "Align the reads", header: true },
		{ name: "Data Analysis", component: Step3, subtitle: "Generate a consensus sequence" },
		{ name: "Data Analysis", component: Step4, subtitle: "Avoid intermediate files" },
		{ name: "The end", component: Conclusion, subtitle: "Conclusion", header: true }
	],
	files: ["primer.bed", "reference.fasta", "reads_R1.fq", "reads_R2.fq"],
	init: `REF_FASTA=reference.fasta
PRIMER_BED=primer.bed`
};
