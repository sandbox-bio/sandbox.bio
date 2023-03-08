// Steps
import Step0 from "./steps/step00.md";
import Step1 from "./steps/step01.md";
import Step2 from "./steps/step02.md";
import Step3 from "./steps/step03.md";
import Step4 from "./steps/step04.md";
import Step5 from "./steps/step05.md";
import Step6 from "./steps/step06.md";
import Step7 from "./steps/step07.md";
import Step8 from "./steps/step08.md";

export const config = {
	id: "ifb-2",
	pwd: "ifb-2",
	listed: false,
	name: "Manipulating files and directories",
	subtitle: `by <a href="https://www.france-bioinformatique.fr/en/home/" target="_blank">French Institute of Bioinformatics</a>`,
	description: "IFB Scenario 2",
	tags: ["unix", "shell", "terminal"],
	tools: ["ls", "date"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Manipulating files and directories", component: Step0 },
		{ name: "Title Step 1", component: Step1 },
		{ name: "Title Step 2", component: Step2 },
		{ name: "Title Step 3", component: Step3 },
		{ name: "Title Step 4", component: Step4 },
		{ name: "Title Step 5", component: Step5 },
		{ name: "Title Step 6", component: Step6 },
		{ name: "Title Step 7", component: Step7 },
		{ name: "Title Step 8", component: Step8 },
	],
	files: [
		"data/ifb-2/absolute_and_relative_paths.png",
	],
};
