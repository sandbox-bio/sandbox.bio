// Steps
import Intro from "./steps/Intro.md";
import Step01 from "./steps/Step01.md";
import Step02 from "./steps/Step02.md";
import Step03 from "./steps/Step03.md";
import Conclusion from "./steps/Conclusion.md";
import Exercise01 from "./exercises/Exercise01.md";

export const config = {
	id: "datepro-01-unix-shell",
	pwd: "datepro-01-unix-shell",
	name: "Data and Text Processing: Unix Shell",
	subtitle: `by <a href="https://webpages.ciencias.ulisboa.pt/~fjcouto/" target="_blank">Francisco M. Couto</a>`,
	description: "Unix Shell",
	tags: ["unix", "shell", "terminal", "script", "file"],
	tools: ["ls", "cat", "tac", "ps", "tr", "mv", "cd", "pwd", "nano", "chmod"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "Unix Shell", component: Step01 },
		{ name: "Data and Script Files", component: Step02 },
		{ name: "Save the Output", component: Step03 },
		{ name: "Conclusion", component: Conclusion },
		{ name: "Exercise", component: Exercise01 }
	],
	files: []
};
