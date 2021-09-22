import Intro from "./steps/Intro.md";
import Exercise from "./steps/Exercise.svelte";

export const config = {
	id: "rosalind",
	terminal: false,
	ide: true,
	listed: false,
	// pwd: "samtools-intro",
	name: "Rosalind Exercises",
	tags: [],
	tools: [],
	difficulty: [],
	steps: [
		{ name: "Rosalind Exercises", component: Intro, subtitle: "Introduction" },
		{ name: "Rosalind Exercises", component: Exercise, subtitle: "Counting DNA Nucleotides", header: true },
		{ name: "Rosalind Exercises", component: Exercise, subtitle: "Transcribing DNA into RNA" }
	]
};
