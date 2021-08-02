// Steps
import Intro from "./steps/Intro.md";
// Data
import Reads1 from "./data/reads_1.fq";
import Reads2 from "./data/reads_2.fq";

export const config = {
	id: "bowtie2-intro",
	name: "Introduction to bowtie2",
	description: "Align DNA sequencing reads from <code>.fastq</code> files to the Lambda phage reference genome.",
	tools: ["bowtie2"],
	difficulty: ["beginner"],
	adapted_from: {
		"name": "bowtie-bio.sourceforge.net/bowtie2",
		"link": "http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example"
	},
	steps: [
		{ name: "Introduction to bowtie2", component: Intro },
	],
	files: [
		{ name: "reads_1.fq", contents: Reads1 },
		{ name: "reads_2.fq", contents: Reads2 },
	],
	// init: async () => {
	// 	console.log(await $CLI.exec("ls /"))
	// }
	//"mv /bowtie2/example/index /shared/data/index"
};
