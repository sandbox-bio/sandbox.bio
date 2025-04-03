// Steps
import Intro from "./steps/Intro.md";
import Conclusion from "./steps/Conclusion.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Step7 from "./steps/Step7.md";

export const config = {
	id: "dna-secrets",
	name: "Variant calling",
	subtitle: `by <a href="https://robert.bio" target="_blank">Robert Aboukhalil</a>`,
	description: "Use variant calling to decode a secret message stored in sequencing data.",
	tags: ["bowtie2", "bcftools"],
	tools: ["samtools", "bcftools", "bowtie2", "echo", "ls", "head"],
	difficulty: ["intermediate"],
	steps: [
		{ name: "Your mission: Variant Calling", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "Data analysis", component: Step2, subtitle: "Align reads to the genome", header: true },
		{ name: "Data analysis", component: Step3, subtitle: "Call variants" },
		{ name: "Data analysis", component: Step4, subtitle: "View the variants" },
		{ name: "Data analysis", component: Step5, subtitle: "Plot twist" },
		{ name: "Data analysis", component: Step6, subtitle: "Extract encoded message" },
		{ name: "Data analysis", component: Step7, subtitle: "Decode the message" },
		{ name: "The end", component: Conclusion, header: true }
	],
	files: ["reads.fq", "morereads.fq"],
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
export REF_FASTA=./reference/lambda_virus.fa;
echo "CGGCGAACAGGCCTAGATTAGGCCCTTCTTCCCGGCGGTG" > secret;`
};
