// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Exercise1 from "./steps/Exercise1.md";
import Exercise2 from "./steps/Exercise2.md";
import Exercise3 from "./steps/Exercise3.md";
import Exercise4 from "./steps/Exercise4.md";
import Exercise5 from "./steps/Exercise5.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "seqkit-intro",
	name: "Wrangle FASTA and FASTQ with SeqKit",
	description: "Explore and wrangle <code>.fasta/.fastq</code> files with SeqKit",
	tags: ["draft", "fasta", "fastq", "seqkit"],
	tools: ["ls", "seqkit"],
	difficulty: ["beginner"],
	new: true,
	steps: [
		{ name: "Wrangle FASTA and FASTQ with SeqKit", component: Intro },
		{ name: "Getting Started", component: Step1, subtitle: "Calculate summary stats", header: true },
		{ name: "Getting Started", component: Step2, subtitle: "Downsample data" },
		{ name: "Getting Started", component: Step3, subtitle: "Extract and filter" },
		{ name: "Getting Started", component: Step4, subtitle: "Extract sequence IDs" },
		{ name: "Getting Started", component: Step5, subtitle: "Update sequence IDs" },
		{ name: "Getting Started", component: Step6, subtitle: "Removing duplicates" },
		{ name: "Exercises", component: Exercise1, subtitle: "Bases per line", header: true },
		{ name: "Exercises", component: Exercise2, subtitle: "Only keep sequence IDs" },
		{ name: "Exercises", component: Exercise3, subtitle: "Searching" },
		{ name: "Exercises", component: Exercise4, subtitle: "Downsample paired-end reads" },
		{ name: "Exercises", component: Exercise5, subtitle: "Replacing FASTA header" },
		{ name: "The End", component: Conclusion, header: true }
	],
	// head -n10000 <(curl -s https://www.mirbase.org/download/hairpin.fa) > hairpins.fa
	// curl -s "https://42basepairs.com/download/r2/genomics-data/reads_NA12878_R1.fastq.gz" | gunzip -c | head -n4000 > NA12878.fastq
	files: ["hairpins.fa", "NA12878.fastq"]
};
