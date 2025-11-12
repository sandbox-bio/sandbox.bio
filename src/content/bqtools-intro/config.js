// Steps
import Intro from "./steps/Intro.md";
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
import Conclusion from "./steps/Conclusion.md";
// Exercises
import Exercise1 from "./exercises/Exercise1.md";
import Exercise2 from "./exercises/Exercise2.md";

export const config = {
	// Metadata
	id: "bqtools-intro",
	name: "Introduction to bqtools",
	subtitle: `by <a href="https://google.com" target="_blank">Noam Teyssier</a>`,
	description: "Exploring bqtools, a toolkit for working with and manipulating binseq files",
	tags: ["bqtools", "binseq", "sequences", "cli"],
	difficulty: ["beginner"],

	// Preload these tools as soon as the page loads
	tools: ["bqtools", "ls", "head", "tail", "sort", "grep", "zcat"],

	// Order of steps
	steps: [
		{ name: "Overview of tutorial", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "CLI help", component: Step2 },
		{ name: "Basic encoding", component: Step3 },
		{ name: "Paired encoding", component: Step4 },
		{ name: "Recursive encoding", component: Step5 },
		{ name: "Decoding", component: Step6 },
		{ name: "Concatenation", component: Step7 },
		{ name: "Grep - Simple", component: Step8 },
		{ name: "Grep - Paired", component: Step9 },
		{ name: "Grep - Pattern Counting", component: Step10 },
		// We use subtitle to define subsections. Click the "1 / 6" button at the bottom to see the effect on the table of contents
		{ name: "Exercises", component: Exercise1, subtitle: "Find non-exons", header: true },
		{ name: "Exercises", component: Exercise2, subtitle: "Find exons in 500kb regions" },
		{ name: "The end", component: Conclusion, header: true }
	],

	// Files within "data/" that you need at runtime.
	// Whenever you update files within "data/", you will need to restart the "./setup.sh" script.
	files: [
		"fastq/sample1_R1.fastq.gz",
		"fastq/sample1_R2.fastq.gz",
		"fastq/sample2_R1.fastq.gz",
		"fastq/sample2_R2.fastq.gz",
		"fastq/sample3_R1.fastq.gz",
		"fastq/sample3_R2.fastq.gz",
		"fastq/sample4_R1.fastq.gz",
		"fastq/sample4_R2.fastq.gz"
	]
};
