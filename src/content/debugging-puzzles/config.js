// Steps
import Intro from "./steps/Intro.md";
import PuzzleTxtSort from "./steps/PuzzleTxtSort.md";
import PuzzleBedSpaces from "./steps/PuzzleBedSpaces.md";
import PuzzleSamQuery from "./steps/PuzzleSamQuery.md";
import PuzzleFastaFilter from "./steps/PuzzleFastaFilter.md";

export const config = {
	id: "debugging-puzzles",
	name: "Debugging Puzzles",
	subtitle: `by <a href="https://robert.bio" target="_blank">Robert Aboukhalil</a>`,
	description: "Debug file format issues that are commonly seen in genomics.",
	tags: ["samtools", "bedtools", "grep"],
	tools: ["head", "grep", "cat", "sort", "samtools", "bedtools", "sed", "vim"],
	difficulty: ["difficult"],
	steps: [
		{ name: "Bioinformatics Debugging Puzzles", component: Intro },
		{ name: "Puzzles", component: PuzzleFastaFilter, subtitle: "The FASTA file that could not be filtered", header: true },
		{ name: "Puzzles", component: PuzzleTxtSort, subtitle: "The TXT file that could not be sorted" },
		{ name: "Puzzles", component: PuzzleSamQuery, subtitle: "The SAM file that could not be queried" },
		{ name: "Puzzles", component: PuzzleBedSpaces, subtitle: "The BED file that could not be merged" }
	],
	// exons.bed is similar to the ones from bedtools-intro but with space instead of tabs
	files: ["exons.bed", "chromosomes.txt", "alignments.bam", "sequences.fa"]
};
