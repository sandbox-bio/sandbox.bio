import Intro from "./steps/Intro.md";

export const config = {
	id: "rosalind",
	terminal: false,
	listed: false,
	// pwd: "samtools-intro",
	name: "Rosalind Exercises",
	description: "Explore, process and manipulate <code>.sam</code> and <code>.bam</code> files with samtools.",
	tags: ["rosalind"],
	tools: [],
	difficulty: [],
	steps: [
		{ name: "Rosalind Exercises", component: Intro, subtitle: "sdf" },
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
	// files: [
	// 	"data/samtools-intro/sample.sam"
	// ],
};
