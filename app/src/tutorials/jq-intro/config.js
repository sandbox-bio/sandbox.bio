// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";

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
		// { name: "The end", component: Conclusion, header: true }
	],
	files: [
		"data/jq-intro/repo.json"
	],
};
