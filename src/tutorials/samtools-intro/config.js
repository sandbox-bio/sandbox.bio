// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";

export const config = {
	id: "samtools-intro",
	name: "Introduction to samtools",
	description: ".",
	tools: ["samtools"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction to samtools", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "Samtools help", component: Step2 },
		{ name: "Samtools utilities", component: Step3, subtitle: "Converting SAM to BAM", header: true },
		{ name: "Samtools utilities", component: Step4, subtitle: "Sort BAM files" },
		{ name: "Samtools utilities", component: Step5, subtitle: "Index BAM files" },
		{ name: "Samtools utilities", component: Step6, subtitle: "Explore BAM files" },
	],
	files: [
		"sample.sam"
	],
};
