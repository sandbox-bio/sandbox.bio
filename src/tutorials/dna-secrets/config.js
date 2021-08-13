// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";

export const config = {
	id: "dna-secrets",
	name: "Decode DNA secrets",
	description: "Use variant calling and genome arithmetic to decode a secret message encoded in sequencing data.",
	tags: ["bowtie2", "bcftools", "bedtools", "samtools"],
	tools: ["samtools/1.10", "bcftools/1.10", "bedtools/2.29.2", "bowtie2/bowtie2-align-s/2.4.2"],
	difficulty: ["intermediate"],
	steps: [
		{ name: "Your mission", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "Steps", component: Step2, subtitle: "Align reads to the genome", header: true },
		{ name: "Steps", component: Intro, subtitle: "Call variants" },
		{ name: "Steps", component: Intro, subtitle: "Plot twist" },
		{ name: "Steps", component: Intro, subtitle: "Extract encoded message" },
		{ name: "Steps", component: Intro, subtitle: "Decode the message" },
	],
	files: [
		"data/dna-secrets/reads.fq"
	],
	init: "REF=/bowtie2/example/index/lambda_virus; REF_FASTA=/bowtie2/example/reference/lambda_virus.fa"
};
