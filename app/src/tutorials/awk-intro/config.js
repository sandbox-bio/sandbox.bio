// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
// import Step5 from "./steps/Step5.md";
// import Step6 from "./steps/Step6.md";
// import Step7 from "./steps/Step7.md";
// import Step8 from "./steps/Step8.md";
// import Step9 from "./steps/Step9.md";
// import Step10 from "./steps/Step10.md";
// import Step11 from "./steps/Step11.md";
// import Step12 from "./steps/Step12.md";
// import Step13 from "./steps/Step13.md";
// import Conclusion from "./steps/Conclusion.md";
// Exercises
import Exercise1 from "./exercises/Exercise1.md";
import Exercise2 from "./exercises/Exercise2.md";
// import Exercise3 from "./exercises/Exercise3.md";
// import Exercise4 from "./exercises/Exercise4.md";
// import Exercise5 from "./exercises/Exercise5.md";

export const config = {
	id: "awk-intro",
	pwd: "awk-intro",
	name: "Intro to awk",
	description: "FIXME:",
	tags: ["awk", "terminal"],
	tools: ["awk"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction to awk", component: Intro },
		{ name: "Filtering data", subtitle: "Extract columns", component: Step1, header: true },
		{ name: "Filtering data", subtitle: "Extract rows", component: Step2 },
		{ name: "Filtering data", subtitle: "Exercise", component: Exercise1 },
		{ name: "Variables", subtitle: "Using variables to count", component: Step3, header: true },
		{ name: "Variables", subtitle: "Auto-initialization", component: Step4 },
		{ name: "Variables", subtitle: "Exercises", component: Exercise2},
		// { name: "The end", component: Conclusion, header: true }
	],
	files: [
		"data/awk-intro/orders.tsv"
	]
}
