// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "bowtie2-intro",
	name: "Sequence alignment with bowtie2",
	icon: "list-nested",
	subtitle: `by <a href="https://langmead-lab.org/" target="_blank">Ben Langmead</a>`,
	description: "Align reads to a reference genome.",
	tags: ["bowtie2", "bcftools"],
	tools: ["bowtie2", "samtools", "bcftools", "ls", "echo", "head"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Sequence alignment with bowtie2", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "Aligning example reads", component: Step2, header: true, subtitle: "Align single-end reads" },
		{ name: "Aligning example reads", component: Step3, subtitle: "Align paired-end reads" },
		{ name: "Aligning example reads", component: Step4, subtitle: "Local alignment" },
		{ name: "Downstream tools", component: Step5, header: true, subtitle: "Use samtools to explore alignments" },
		{ name: "Downstream tools", component: Step6, subtitle: "Use bcftools to call variants" },
		{ name: "The end", component: Conclusion, header: true }
	],
	files: ["reads_1.fq", "reads_2.fq", "longreads.fq"],
	assets: [
		"reference/lambda_virus.fa",
		"index/lambda_virus.1.bt2",
		"index/lambda_virus.2.bt2",
		"index/lambda_virus.3.bt2",
		"index/lambda_virus.4.bt2",
		"index/lambda_virus.rev.1.bt2",
		"index/lambda_virus.rev.2.bt2"
	],
	init: `export REF=./index/lambda_virus;
export REF_FASTA=./reference/lambda_virus.fa;`
};
