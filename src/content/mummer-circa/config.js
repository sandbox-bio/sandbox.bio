// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Conclusion from "./steps/Conclusion.md";
// Exercises
// import Exercise1 from "./exercises/Exercise1.md";
// import Exercise2 from "./exercises/Exercise2.md";

export const config = {
	// Metadata
	id: "mummer-circa",
	name: "Aligning genomes with MUMmer",
	subtitle: `by <a href="https://www.linkedin.com/in/marianattestad/" target="_blank">Maria Nattestad</a>`,
	description: "Aligning genomes with MUMmer",
	tags: ["mummer", "circa"],
	difficulty: ["intermediate"],

	// Preload these tools as soon as the page loads
	tools: ["nucmer", "show-coords", "ls", "echo", "awk", "samtools", "cat"],

	// Order of steps
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "Bedtools help", component: Step2 },
		// We use subtitle to define subsections. Click the "1 / 6" button at the bottom to see the effect on the table of contents
		// { name: "Exercises", component: Exercise1, subtitle: "Find non-exons", header: true },
		// { name: "Exercises", component: Exercise2, subtitle: "Find exons in 500kb regions" },
		{ name: "The end", component: Conclusion, header: true }
	],

	// Files within "data/" that you need at runtime.
	// Whenever you update files within "data/", you will need to restart the "./setup.sh" script.
	files: ["Podan2_AssemblyScaffoldsmt.fa", "CBS415.72m.nice_mt.fa"]
};
