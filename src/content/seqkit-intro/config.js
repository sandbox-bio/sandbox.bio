// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Exercise1 from "./steps/Exercise1.md";
import Exercise2 from "./steps/Exercise2.md";
import Exercise3 from "./steps/Exercise3.md";
import Exercise4 from "./steps/Exercise4.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "seqkit-intro",
	name: "Wrangle FASTA and FASTQ with SeqKit",
	icon: "wrench-adjustable",
	description: "Explore and wrangle <code>.fasta</code> and <code>.fastq</code> files with SeqKit",
	tags: ["fasta", "fastq", "seqkit"],
	tools: ["ls", "seqkit"],
	difficulty: ["beginner"],
	new: true,
	steps: [
		{ name: "Wrangle FASTA and FASTQ with SeqKit", component: Intro },
		{ name: "Getting Started", component: Step1, subtitle: "Calculate summary stats", header: true },
		{ name: "Getting Started", component: Step2, subtitle: "Downsample data" },
		{ name: "Getting Started", component: Step3, subtitle: "Extract and filter" },
		{ name: "Getting Started", component: Step4, subtitle: "Search FASTA/FASTQ files" },
		{ name: "Getting Started", component: Step5, subtitle: "Remove duplicate sequences" },
		{ name: "Exercises", component: Exercise1, subtitle: "Format sequences", header: true },
		{ name: "Exercises", component: Exercise2, subtitle: "Only keep sequence IDs" },
		{ name: "Exercises", component: Exercise3, subtitle: "Search for sequences" },
		{ name: "Exercises", component: Exercise4, subtitle: "Renaming sequence names" },
		{ name: "The End", component: Conclusion, header: true }
	],
	// head -n10000 <(curl -s https://www.mirbase.org/download/hairpin.fa) > hairpins.fa
	// head -n9 <(curl -s https://www.mirbase.org/download/hairpin.fa) >> hairpins.fa  # adding duplicates for Step 5
	// curl -s "https://42basepairs.com/download/r2/genomics-data/reads_NA12878_R1.fastq.gz" | gunzip -c | head -n4000 > NA12878.fastq
	files: ["hairpins.fa", "NA12878.fastq", "ids.txt"]
};
