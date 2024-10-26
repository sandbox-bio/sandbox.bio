// Steps
import Intro from "./steps/Intro.md";
import Step01 from "./steps/Step01.md";
import Step02 from "./steps/Step02.md";
import Step03 from "./steps/Step03.md";
import Step04 from "./steps/Step04.md";
import Step05 from "./steps/Step05.md";
import Step06 from "./steps/Step06.md";
import Step07 from "./steps/Step07.md";
import Conclusion from "./steps/Conclusion.md";
import Exercise01 from "./exercises/Exercise01.md";

export const config = {
	id: "datepro-07-text-processing",
	pwd: "datepro-07-text-processing",
	name: "Chapter 4: Text Processing",
	subtitle: `by <a href="https://webpages.ciencias.ulisboa.pt/~fjcouto/" target="_blank">Francisco M. Couto</a>`,
	description: "Pattern matching and regular expressions to extract information from text.",
	tags: ["Pattern matching", "Regular expressions", "Tokenization", "Entity recognition", "Relation extraction"],
	tools: ["grep", "diff", "sed", "echo", "sort", "tr", "wc"],
	difficulty: ["advanced"],
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "Pattern Matching", component: Step01 },
		{ name: "Regular Expressions", component: Step02 },
		{ name: "Position", component: Step03 },
		{ name: "Tokenization", component: Step04 },
		{ name: "Entity recognition", component: Step05 },
		{ name: "Pattern File", component: Step06 },
		{ name: "Relation Extraction", component: Step07 },
		{ name: "Conclusion", component: Conclusion },
		{ name: "Exercise", component: Exercise01 }
	],
	files: ["chebi_27732.txt"]
};
