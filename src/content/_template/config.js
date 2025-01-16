// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Conclusion from "./steps/Conclusion.md";
// Exercises
import Exercise1 from "./exercises/Exercise1.md";
import Exercise2 from "./exercises/Exercise2.md";

export const config = {
	// Metadata
	id: "_template",
	name: "An example tutorial",
	subtitle: `by <a href="https://google.com" target="_blank">Your Name</a>`,
	description: "Do things and learn things.",
	tags: ["bedtools", "igv"],
	difficulty: ["beginner"],

	// Preload these tools as soon as the page loads
	tools: ["bedtools", "ls", "head", "tail", "sort"],

	// Order of steps
	steps: [
		{ name: "An example tutorial", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "Bedtools help", component: Step2 },
		// We use subtitle to define subsections. Click the "1 / 6" button at the bottom to see the effect on the table of contents
		{ name: "Exercises", component: Exercise1, subtitle: "Find non-exons", header: true },
		{ name: "Exercises", component: Exercise2, subtitle: "Find exons in 500kb regions" },
		{ name: "The end", component: Conclusion, header: true }
	],

	// Files within "data/" that you need at runtime.
	// Whenever you update files within "data/", you will need to restart the "./setup.sh" script.
	files: ["cpg.bed", "exons.bed", "fHeart-DS15839.bed", "fHeart-DS16621.bed", "fSkin-DS19745.bed", "gwas.bed", "hesc.chromHmm.bed", "genome.txt"]
};
