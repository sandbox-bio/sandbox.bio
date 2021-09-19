import Intro from "./steps/Intro.md";
import Exercise from "./steps/Exercise.md";

export const config = {
	id: "rosalind",
	terminal: false,
	ide: true,
	listed: false,
	// pwd: "samtools-intro",
	name: "Rosalind Exercises",
	description: "Explore, process and manipulate <code>.sam</code> and <code>.bam</code> files with samtools.",
	tags: [],
	tools: [],
	difficulty: [],
	steps: [
		{ name: "Rosalind Exercises", component: Intro, subtitle: "Welcome" },
		{ name: "Rosalind Exercises", component: Exercise, subtitle: "Counting DNA Nucleotides", header: true },
		{ name: "Rosalind Exercises", component: Exercise, subtitle: "Transcribing DNA into RNA" }
	]
};
