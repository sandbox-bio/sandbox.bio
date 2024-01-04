import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";

export const config = {
	id: "blast-intro",
	name: "Sequence alignment with BLAST",
	subtitle: `by <a href="https://shawntoneil.com/" target="_blank">Shawn T. O'Neil</a>`,
	description: "Use BLAST to align DNA and protein sequences",
	tags: ["blastn", "blastp"],
	tools: ["makeblastdb", "blastdbcmd", "blastp", "blastn"],
	difficulty: ["beginner"],

	steps: [
		{ name: "Sequence alignment with BLAST", component: Intro },
		{ name: "BLAST Types", component: Step1 },
		{ name: "BLAST Databases", component: Step2 },
		{ name: "Run BLAST", subtitle: "Create a BLAST database", component: Step3, header: true },
		{ name: "Run BLAST", subtitle: "Run blastp", component: Step4 },
		{ name: "Run BLAST", subtitle: "Interpret the results", component: Step5 },
		{ name: "Exercises", subtitle: "Run BLAST", component: Step6 }
	],
	files: ["orf_trans.fasta", "p450s.fasta"]
};
