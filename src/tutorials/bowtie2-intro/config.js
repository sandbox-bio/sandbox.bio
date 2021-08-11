// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "bowtie2-intro",
	name: "Introduction to bowtie2",
	description: "Align DNA sequencing reads from <code>.fastq</code> files to a reference genome",
	tags: ["bowtie2", "samtools"],
	tools: ["samtools/1.10", "bowtie2/bowtie2-align-s/2.4.2"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction to bowtie2", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "Aligning example reads", component: Step2, header: true, subtitle: "Align single-end reads" },
		{ name: "Aligning example reads", component: Step3, subtitle: "Align paired-end reads" },
		{ name: "Aligning example reads", component: Step4, subtitle: "Local alignment" },
		{ name: "Using samtools downstream", component: Step5, subtitle: "Using samtools" },
		{ name: "The end", component: Conclusion, header: true }
	],
	files: [
		"data/bowtie2-intro/reads_1.fq",
		"data/bowtie2-intro/reads_2.fq",
		"data/bowtie2-intro/longreads.fq",
		"data/bowtie2-intro/r1.fq",
	],
	init: "REF=/bowtie2/example/index/lambda_virus; REF_FA=/bowtie2/example/reference/lambda_virus.fa"
};
