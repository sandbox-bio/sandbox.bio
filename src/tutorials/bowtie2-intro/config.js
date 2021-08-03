// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";

export const config = {
	id: "bowtie2-intro",
	name: "Introduction to bowtie2",
	description: "Align DNA sequencing reads from <code>.fastq</code> files to the Lambda phage reference genome.",
	tools: ["bowtie2"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction to bowtie2", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "Aligning example reads", component: Step2, header: true, subtitle: "Align single-end reads" },
		{ name: "Aligning example reads", component: Step3, subtitle: "Align paired-end reads" },
		{ name: "Aligning example reads", component: Step4, subtitle: "Local alignment" },
		{ name: "Using samtools downstream", component: Step5, subtitle: "Using samtools" },
	],
	files: [
		"reads_1.fq",
		"reads_2.fq",
		"longreads.fq"
	],
	init: "REF=/bowtie2/example/index/lambda_virus"
};
