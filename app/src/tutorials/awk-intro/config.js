// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step4_2 from "./steps/Step4_2.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Conclusion from "./steps/Conclusion.md";
// Exercises
import Exercise1 from "./exercises/Exercise1.md";
import Exercise2 from "./exercises/Exercise2.md";
import Exercise3 from "./exercises/Exercise3.md";

export const config = {
	id: "awk-intro",
	pwd: "awk-intro",
	name: "Data exploration with awk",
	subtitle: `by <a href="https://robert.bio" target="_blank">Robert Aboukhalil</a>`,
	description: "Filter, extract and transform data files using <code>awk</code>.",
	tags: ["awk", "terminal"],
	tools: ["awk"],
	difficulty: ["intermediate"],
	steps: [
		{ name: "Data exploration with awk", component: Intro },
		{ name: "Filtering data", subtitle: "Extract columns", component: Step1, header: true },
		{ name: "Filtering data", subtitle: "Extract rows", component: Step2 },
		{ name: "Filtering data", subtitle: "Exercise", component: Exercise1 },
		{ name: "Variables", subtitle: "Using variables to count sums", component: Step3, header: true },
		{ name: "Variables", subtitle: "Auto-initialization of variables", component: Step4 },
		{ name: "Variables", subtitle: "Passing variables from the outside", component: Step4_2 },
		{ name: "Variables", subtitle: "Exercise", component: Exercise2 },
		{ name: "Arrays", subtitle: "Using arrays to tally burritos", component: Step5, header: true },
		{ name: "Arrays", subtitle: "Using arrays to tally costs", component: Step6 },
		{ name: "Arrays", subtitle: "Exercise", component: Exercise3 },
		{ name: "The end", component: Conclusion, header: true }
	],
	files: [
		"data/awk-intro/orders.tsv"
	]
}
