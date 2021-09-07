// Steps
import Intro from "./steps/Intro.md";

export const config = {
	id: "jq-intro",
	name: "Introduction to jq",
	description: "Filter and wrangle <code>JSON</code> with jq.",
	tags: ["jq", "terminal"],
	tools: ["jq/1.6"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction to jq", component: Intro },
		// { name: "The data", component: Step1 },
		// { name: "Samtools help", component: Step2 },
		// { name: "Samtools utilities", component: Step3, subtitle: "Converting SAM to BAM", header: true },
		// { name: "Samtools utilities", component: Step4, subtitle: "Sort BAM files" },
		// { name: "Samtools utilities", component: Step5, subtitle: "Index BAM files" },
		// { name: "Explore BAM files", component: Step6, subtitle: "Scrutinize alignments", header: true },
		// { name: "Explore BAM files", component: Step7, subtitle: "Inspect the header" },
		// { name: "Explore BAM files", component: Step8, subtitle: "Capture the flag" },
		// { name: "The end", component: Conclusion, header: true }
	],
	files: [
		// "data/samtools-intro/sample.sam"
	],
};
