// Steps
import Intro from "./steps/Intro.md";

export const config = {
	id: "hidden-message",
	name: "Decode data stored in DNA",
	description: "TODO",
	tools: ["bowtie2", "bcftools", "IGV"],
	difficulty: ["intermediate"],
	steps: [
		{ name: "Your mission", component: Intro },
		{ name: "The data", component: Intro },
		{ name: "Steps", component: Intro, subtitle: "Align reads to the genome", header: true },
		{ name: "Steps", component: Intro, subtitle: "Call variants" },
		{ name: "Steps", component: Intro, subtitle: "Extract encoded message" },
		{ name: "Steps", component: Intro, subtitle: "Decode the message" },
	],
	files: [
		"data/hidden-message/reads.fq"
	],
	init: "REF=/bowtie2/example/index/lambda_virus; REF_FASTA=/bowtie2/example/reference/lambda_virus.fa"
};
