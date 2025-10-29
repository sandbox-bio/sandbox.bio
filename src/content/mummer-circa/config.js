// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";

export const config = {
	// Metadata
	id: "mummer-circa",
	icon: "bezier2",
	date: "October 2025",
	name: "Aligning genomes with MUMmer",
	subtitle: `by <a href="https://www.linkedin.com/in/marianattestad/" target="_blank">Maria Nattestad</a>`,
	description: "Use MUMmer to align two bacterial genomes and visualize the results with Circa.",
	tags: ["mummer", "circa"],
	difficulty: ["intermediate"],

	// Preload these tools
	tools: ["nucmer", "show-coords", "samtools", "awk", "cat", "echo", "head", "less"],

	// Order of steps
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "Genome Alignment", component: Step1 },
		{ name: "Data Processing", component: Step2 },
		{ name: "Reference Files", component: Step3 },
		{ name: "Visualization", component: Step4 }
	],

	// Files within "data/" that the tutorial needs at runtime.
	// Whenever you update files within "data/", you will need to restart the "./setup.sh" script.
	files: ["H_pylori26695_Eslice.fasta", "H_pyloriJ99_Eslice.fasta"]
};
