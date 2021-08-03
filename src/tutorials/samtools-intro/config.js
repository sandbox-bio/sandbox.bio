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
	id: "samtools-intro",
	name: "Introduction to samtools",
	description: "Explore, process and manipulate <code>.sam</code> and <code>.bam</code> files with samtools",
	tools: ["samtools"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction to samtools", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "Samtools help", component: Step2 },
		{ name: "Samtools utilities", component: Step3, subtitle: "Converting SAM to BAM", header: true },
		{ name: "Samtools utilities", component: Step4, subtitle: "Sort BAM files" },
		{ name: "Samtools utilities", component: Step5, subtitle: "Index BAM files" },
		{ name: "Explore BAM files", component: Step6, subtitle: "Scrutinize alignments", header: true },
		{ name: "Explore BAM files", component: Step7, subtitle: "Inspect the header" },
		{ name: "Explore BAM files", component: Step8, subtitle: "Capture the flag" },
		{ name: "The end", component: Conclusion, header: true }
	],
	files: [
		"sample.sam"
	],
};
