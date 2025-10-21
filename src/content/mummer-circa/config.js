// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";

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

	// Preload these tools
	tools: ["nucmer", "show-coords", "samtools", "awk", "cat", "echo"],

	// Order of steps
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "Genome Alignment", component: Step1 },
		{ name: "Data Processing", component: Step2 },
		{ name: "Reference Files", component: Step3 },
		{ name: "Visualization", component: Step4 }
	],

	// Files within "data/" that you need at runtime.
	// Whenever you update files within "data/", you will need to restart the "./setup.sh" script.
	files: ["H_pylori26695_Eslice.fasta", "H_pyloriJ99_Eslice.fasta"]
};
