// Steps
import Intro from "./steps/Intro.md";
import Step01 from "./steps/Step01.md";
import Step02 from "./steps/Step02.md";
import Step03 from "./steps/Step03.md";
import Step04 from "./steps/Step04.md";
import Step05 from "./steps/Step05.md";
import Step06 from "./steps/Step06.md";
import Step07 from "./steps/Step07.md";
import Step08 from "./steps/Step08.md";
import Step09 from "./steps/Step09.md";
import Conclusion from "./steps/Conclusion.md";
import Exercise01 from "./exercises/Exercise01.md";

export const config = {
	id: "datepro-08-semantic-processing",
	pwd: "datepro-08-semantic-processing",
	name: "Chapter 5: Semantic Processing",
	subtitle: `by <a href="https://webpages.ciencias.ulisboa.pt/~fjcouto/" target="_blank">Francisco M. Couto</a>`,
	description: "Exploring biomedical ontologies to identify relevant entities in text.",
	tags: ["Ontologies", "Classes", "Ancestors", "Lexicon", "Entity Linking"],
	tools: ["xmllint", "curl",  "xargs", "gunzip", "unzip", "echo", "while", "sort", "tr", "wc"],
	difficulty: ["advanced"],
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "Classes", component: Step01 },
		{ name: "URIs and Labels", component: Step02 },
		{ name: "Synonyms", component: Step03 },
		{ name: "Parent Classes", component: Step04 },
		{ name: "Ancestors", component: Step05 },
		{ name: "My Lexicon", component: Step06 },
		{ name: "Generic Lexicon", component: Step07 },
		{ name: "Entity Linking", component: Step08 },
		{ name: "Large lexicons", component: Step09 },
		{ name: "Conclusion", component: Conclusion },
		{ name: "Exercise", component: Exercise01 }
	],
	files: [
		"chebi_lite.owl.gz",
		"doid.owl.gz"
        	]
};
