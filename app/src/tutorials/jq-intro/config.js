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

export const config = {
	id: "jq-intro",
	name: "Introduction to jq",
	description: "Filter and wrangle <code>JSON</code> with jq.",
	tags: ["jq", "terminal"],
	tools: ["jq/1.6"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction to jq", component: Intro },
		{ name: "Select Elements", component: Step1 },
		{ name: "Select Arrays", component: Step2 },
		{ name: "Putting Elements in an Array", component: Step3, subtitle: "Converting SAM to BAM", header: true },
		{ name: "Select Multiple Fields", component: Step4, subtitle: "Sort BAM files" },
		{ name: "Putting Elements Into an Object", component: Step5, subtitle: "Index BAM files" },
		{ name: "Sorting and Counting", component: Step6, subtitle: "Scrutinize alignments", header: true },
		{ name: "Maps and Selects", component: Step7, subtitle: "Inspect the header" },
		{ name: "In Review", component: Step8, subtitle: "Capture the flag" },
		// { name: "The end", component: Conclusion, header: true }
	],
	files: [
		"data/jq-intro/repo.json",
		"data/jq-intro/issues.json",
		"data/jq-intro/issue.json"
	],
};
