// Steps
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Step7 from "./steps/Step7.md";
import Step8 from "./steps/Step8.md";
import Step9 from "./steps/Step9.md";
import Step10 from "./steps/Step10.md";
import Step11 from "./steps/Step11.md";
import Step12 from "./steps/Step12.md";
// import Exercise1 from "./steps/Exercise1.md";
// import Exercise2 from "./steps/Exercise2.md";
// import Exercise3 from "./steps/Exercise3.md";
// import Exercise4 from "./steps/Exercise4.md";
// import Exercise5 from "./steps/Exercise5.md";
// import Exercise6 from "./steps/Exercise6.md";
// import Exercise7 from "./steps/Exercise7.md";
// import Exercise8 from "./steps/Exercise8.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "fqgrep-intro",
	name: "Pattern searching in FASTQ files with fqgrep",
	subtitle: `by <a href="https://github.com/nh13" target="_blank">Nils Homer</a>`,
	description: "Learn to search for sequence patterns in FASTQ files using fqgrep, a fast grep alternative designed for sequencing data.",
	tags: ["FASTQ", "grep", "pattern matching", "paired-end", "QC"],
	tools: ["fqgrep"],
	difficulty: ["beginner", "intermediate"],
	steps: [
		{ name: "Introduction", component: Step1 },
		{ name: "Basic Pattern Matching", component: Step2 },
		{ name: "Color and Progress", component: Step3 },
		{ name: "Multiple Patterns", component: Step4 },
		{ name: "Inverted Matching", component: Step5 },
		{ name: "Paired-End Reads", component: Step6 },
		{ name: "Reverse Complement", component: Step7 },
		{ name: "Regular Expressions", component: Step8 },
		{ name: "IUPAC Ambiguity Codes", component: Step9 },
		{ name: "Read Name Filtering", component: Step10 },
		{ name: "Protein Sequences", component: Step11 },
		{ name: "Performance", component: Step12 }
		// { name: "Adapter Check", component: Exercise1 },
		// { name: "Tn5 Detection", component: Exercise2 },
		// { name: "Degenerate Barcodes", component: Exercise3 },
		// { name: "Paired-End Screening", component: Exercise4 },
		// { name: "Extract Reads by Name", component: Exercise5 },
		// { name: "Protein Motif Search", component: Exercise6 },
		// { name: "IUPAC Search", component: Exercise7 },
		// { name: "Thorough Contaminant Screen", component: Exercise8 },
		{ name: "Conclusion", component: Conclusion },
	],
	files: ["reads.fastq", "reads_R1.fastq", "reads_R2.fastq", "adapters.txt", "barcodes.txt", "contaminants.txt", "read_names.txt", "proteins.fastq"]
};
