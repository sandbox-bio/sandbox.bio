// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "terminal-basics",
	pwd: "terminal-basics",
	name: "Terminal Basics",
	subtitle: `by <a href="https://robert.bio" target="_blank">Robert Aboukhalil</a>`,
	description: "Get up to speed with the terminal. Start here if you're new to bioinformatics.",
	tags: ["terminal"],
	tools: ["echo", "ls", "head", "tail", "grep", "wc", "cat"],
	difficulty: ["beginner"],
	divider: "Terminal",
	steps: [
		{ name: "Terminal Basics", component: Intro },
		{ name: "The Basics", component: Step1, subtitle: "Navigate the terminal", header: true },
		{ name: "The Basics", component: Step2, subtitle: "Preview files" },
		{ name: "The Basics", component: Step3, subtitle: "Filter files" },
		{ name: "Pipelines", component: Step4, subtitle: "Pipes", header: true },
		{ name: "Pipelines", component: Step5, subtitle: "Output to a file" },
		{ name: "Pipelines", component: Step6, subtitle: "Environment variables" },
		{ name: "The end", component: Conclusion, header: true }
	],
	files: [
		"data/terminal-basics/orders.tsv",
		"data/terminal-basics/ref.fa",
		"data/terminal-basics/ref.fa.bak"
	]
};
