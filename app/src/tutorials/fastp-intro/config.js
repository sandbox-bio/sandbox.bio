// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Step5 from "./steps/Step5.md";
import Step6 from "./steps/Step6.md";
import Step7 from "./steps/Step7.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "fastp-intro",
	pwd: "fastp-intro",
	name: "DNA sequencing QC with fastp",
	description: "Evaluate the quality of a sequencing run by running <code>fastp</code> on your FASTQ files.",
	tags: ["fastp", "QC", "sequencing"],
	tools: ["fastp", "jq"],
	difficulty: ["beginner"],
	steps: [
		{ name: "DNA sequencing QC", component: Intro },
		{ name: "The data", component: Step1 },
		{ name: "QC Reports", subtitle: "Basic reports", component: Step2, header: true },
		{ name: "QC Reports", subtitle: "HTML reports", component: Step3 },
		{ name: "QC Reports", subtitle: "JSON reports", component: Step4 },
		{ name: "Filter Options", subtitle: "Output filtered data", component: Step5, header: true },
		{ name: "Filter Options", subtitle: "Filter out short reads", component: Step6 },
		{ name: "Filter Options", subtitle: "Trim low quality bases", component: Step7 },
		{ name: "The end", component: Conclusion, header: true }
	],
	files: [
		"data/fastp-intro/HG004_R1.fastq.gz",
		"data/fastp-intro/HG004_R2.fastq.gz"
	]
}
