// Steps
import Intro from "./steps/Intro.md";
import PuzzleCsvSort from "./steps/PuzzleCsvSort.md";
import PuzzleBedSpaces from "./steps/PuzzleBedSpaces.md";
import PuzzleSamQuery from "./steps/PuzzleSamQuery.md";

export const config = {
	id: "debugging-puzzles",
	pwd: "debugging-puzzles",
	name: "Bioinformatics Debugging Puzzles",
	subtitle: `by <a href="https://robert.bio" target="_blank">Robert Aboukhalil</a>`,
	description: "TODO",
	tags: ["TODO"],
	tools: [],
	difficulty: ["difficult"],
	steps: [
		{ name: "Bioinformatics Debugging Puzzles", component: Intro },
		{ name: "Puzzles", component: PuzzleCsvSort, subtitle: "The CSV file that could not be sorted", header: true },
		{ name: "Puzzles", component: PuzzleBedSpaces, subtitle: "The BED file that could not be merged" },
		{ name: "Puzzles", component: PuzzleSamQuery, subtitle: "The SAM file that could not be queried" },
		// { name: "The end", component: Conclusion, subtitle: "In Review", header: true }
	],
	files: [
		"data/debugging-puzzles/exons.bed", // similar to the ones from bedtools-intro but with space instead of tabs
		"data/debugging-puzzles/chromosomes.csv",
		"data/debugging-puzzles/alignments.bam",
	],
};
