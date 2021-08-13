// Steps
import Intro from "./steps/Intro.md";

export const config = {
	id: "hidden-message",
	name: "Decode DNA secrets",
	description: "Use variant calling and genome arithmetic to decode a secret message encoded in sequencing data.",
	tags: ["bowtie2", "bcftools", "bedtools", "samtools"],
	tools: ["samtools/1.10", "bcftools/1.10", "bedtools/2.29.2", "bowtie2/bowtie2-align-s/2.4.2"],
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
