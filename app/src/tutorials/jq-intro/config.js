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

export const config = {
	id: "jq-intro",
	name: "Introduction to jq",
	description: "Filter and wrangle <code>JSON</code> with jq.",
	tags: ["jq", "terminal"],
	tools: ["jq/1.6"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction to jq", component: Intro },
		{ name: "Filter data", component: Step1, subtitle: "Select Elements", header: true },
		{ name: "Filter data", component: Step2, subtitle: "Select Arrays" },
		{ name: "Filter data", component: Step3, subtitle: "Putting Elements in an Array" },
		{ name: "Filter data", component: Step4, subtitle: "Select Multiple Fields" },
		{ name: "Filter data", component: Step5, subtitle: "Putting Elements Into an Object" },
		{ name: "Summarize data", component: Step6, subtitle: "Sorting and Counting", header: true },
		{ name: "Summarize data", component: Step7, subtitle: "Pipes and Filters" },
		{ name: "Summarize data", component: Step8, subtitle: "Maps and Selects" },
		{ name: "The end", component: Step9, subtitle: "In Review", header: true }
	],
	files: [
		"data/jq-intro/repo.json",
		"data/jq-intro/issues.json",
		"data/jq-intro/issue.json"
	],
};
