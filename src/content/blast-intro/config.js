import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Step7 from "./steps/Step7.md";
import Step8 from "./steps/Step8.md";

export const config = {
	id: "blast-intro",
	name: "Sequence alignment with BLAST",
	description: "Use BLAST to align DNA and protein sequences.",
	tags: ["blastn", "blastp"],
	tools: ["makeblastdb", "blastp", "blastn", "blast_formatter", "blastdbcmd", "wc"],
	difficulty: ["beginner"],
	new: true,
	steps: [
		{ name: "Sequence alignment with BLAST", component: Intro },
		{ name: "BLAST Types", component: Step1 },
		{ name: "BLAST Databases", component: Step2 },
		{ name: "Run BLAST", subtitle: "Create a BLAST database", component: Step3, header: true },
		{ name: "Run BLAST", subtitle: "Run blastp", component: Step4 },
		{ name: "Run BLAST", subtitle: "Interpret the results", component: Step5 },
		{ name: "Exercises", subtitle: "Run BLAST", component: Step6, header: true },
		{ name: "Exercises", subtitle: "Convert formats", component: Step7 },
		{ name: "Exercises", subtitle: "Extract FASTA records", component: Step8 }
	],
	files: ["orf_trans.fasta", "p450s.fasta", "ids.txt"]
};
