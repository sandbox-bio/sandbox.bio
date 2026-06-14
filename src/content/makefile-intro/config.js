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
	id: "makefile-intro",
	name: "Intro to Makefiles",
	subtitle: `by <a href="https://github.com/muppetjones" target="_blank">Stephen Bush</a>`,
	description: "Learn the basics of makefile automation.",
	tags: ["terminal", "shell", "make", "makefile", "vim"],
	difficulty: ["beginner"],

	// Preload these tools as soon as the page loads
	tools: ["make", "ls", "head", "tail", "sort", "vim", "cat", "less"],

	// Order of steps
	steps: [
		{ name: "An Intro to Makefiles", component: Intro },
		{ name: "A Short VIM Refresher", component: Step1 },
		{ name: "Target Practice", component: Step2 },
		{ name: "Fake it, then make it", component: Step3 },
		// The bash-script tutoraial has a text editor, but I don't see that tutorial
		// in the git repo, so for now, introduce the user to vim.
		// { name: "Exercises", component: Exercise1, subtitle: "Find non-exons", header: true },
		// { name: "Exercises", component: Exercise2, subtitle: "Find exons in 500kb regions" },
		{ name: "The end", component: Conclusion, header: true }
	],

	// Files within "data/" that you need at runtime.
	// Whenever you update files within "data/", you will need to restart the "./setup.sh" script.
	files: [
		"Makefile",
	]
};
