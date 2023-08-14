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
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "jq-intro",
	name: "JSON wrangling with jq",
	subtitle: `by <a href="https://adamgordonbell.com/" target="_blank">Adam Gordon Bell</a>`,
	description: "Filter and wrangle JSON files on the command-line using <code>jq</code>.",
	tags: ["jq", "terminal"],
	tools: ["jq", "echo"],
	difficulty: ["beginner"],
	steps: [
		{ name: "JSON wrangling with jq", component: Intro },
		{ name: "Filter data", component: Step1, subtitle: "Select Elements", header: true },
		{ name: "Filter data", component: Step2, subtitle: "Select Arrays" },
		{ name: "Filter data", component: Step3, subtitle: "Putting Elements in an Array" },
		{ name: "Filter data", component: Step4, subtitle: "Select Multiple Fields" },
		{ name: "Filter data", component: Step5, subtitle: "Putting Elements Into an Object" },
		{ name: "Summarize data", component: Step6, subtitle: "Sorting and Counting", header: true },
		{ name: "Summarize data", component: Step7, subtitle: "Pipes and Filters" },
		{ name: "Summarize data", component: Step8, subtitle: "Maps and Selects" },
		{ name: "The end", component: Conclusion, subtitle: "In Review", header: true }
	],
	files: ["data/jq-intro/repo.json", "data/jq-intro/issues.json", "data/jq-intro/issue.json"]
};
