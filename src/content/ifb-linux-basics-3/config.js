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
import Step9 from "./steps/step09.md";

export const config = {
	id: "ifb-linux-basics-3",
	pwd: "ifb-linux-basics-3",
	listed: false,
	name: "Manipulating files and directories",
	subtitle: `by <a href="https://www.france-bioinformatique.fr/en/home/" target="_blank">French Institute of Bioinformatics</a>`,
	description: "IFB Scenario 3",
	tags: ["unix", "shell", "terminal"],
	tools: ["ls", "date"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction", component: Step0 },
		{ name: "Manipulating data", component: Step1 },
		{ name: "Digging through a large file", component: Step2 },
		{ name: "The begining of a file", component: Step3 },
		{ name: "The end of a file", component: Step4 },
		{ name: "Counting words and lines in a file", component: Step5 },
		{ name: "Searching patterns", component: Step6 },
		{ name: "Extracting colums", component: Step7 },
		{ name: "Some useful commands and tips", component: Step8 },
		{ name: "Congratulations", component: Step9 }
	],
	files: ["MACS2.csv", "NC_009089.bed", "NC_009089.fasta", "SAOUHSC.fasta", "SAOUHSC.bed"]
};
