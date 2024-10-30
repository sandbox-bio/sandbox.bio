// Steps
import Intro from "./steps/Intro.md";
import Step01 from "./steps/Step01.md";
import Step02 from "./steps/Step02.md";
import Step03 from "./steps/Step03.md";
import Conclusion from "./steps/Conclusion.md";
import Exercise01 from "./exercises/Exercise01.md";

export const config = {
	id: "datepro-04-task-repetition",
	pwd: "datepro-04-task-repetition",
	name: "Section 3.6: Task Repetition",
	subtitle: `by <a href="https://webpages.ciencias.ulisboa.pt/~fjcouto/" target="_blank">Francisco M. Couto</a>`,
	description: "Streamlining UniProt data retrieval by automating the download of XML files for multiple protein lists in a single step",
	tags: ["proteins", "download", "lists of arguments", "XML", "UniProt"],
	tools: ["xargs", "curl"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "Assembly line", component: Step01 },
		{ name: "Download files", component: Step02 },
		{ name: "Create Script", component: Step03 },
		{ name: "Conclusion", component: Conclusion },
		{ name: "Exercise", component: Exercise01 }
	],
	files: ["chebi_27732_xrefs_UniProt_relevant_identifiers.csv"]
};
