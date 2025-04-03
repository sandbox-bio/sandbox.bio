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
	name: "BAM parsing with samtools",
	subtitle: `by <a href="http://quinlanlab.org/" target="_blank">Aaron Quinlan</a>`,
	description: "Explore and wrangle <code>.sam/.bam</code> files with samtools.",
	tags: ["samtools", "IGV"],
	tools: ["samtools", "ls", "head"],
	difficulty: ["beginner"],
	steps: [
		{ name: "BAM parsing with samtools", component: Intro },
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
	files: ["sample.sam"]
};
