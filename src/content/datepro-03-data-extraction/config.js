// Steps
import Intro from "./steps/Intro.md";
import Step01 from "./steps/Step01.md";
import Step02 from "./steps/Step02.md";
import Step03 from "./steps/Step03.md";
import Conclusion from "./steps/Conclusion.md";
import Exercise01 from "./exercises/Exercise01.md";

export const config = {
	id: "datepro-03-data-extraction",
	pwd: "datepro-03-data-extraction",
	name: "Section 3.5: Data Extraction",
	subtitle: `by <a href="https://webpages.ciencias.ulisboa.pt/~fjcouto/" target="_blank">Francisco M. Couto</a>`,
	description: "Filtering and extracting key data from CSV files, with a focus on proteins linked to potential disease associations",
	tags: ["filter", "extract", "pattern", "CSV"],
	tools: ["grep", "cut"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "Single pattern", component: Step01 },
		{ name: "Multiple patterns", component: Step02 },
		{ name: "Extract column", component: Step03 },
		{ name: "Conclusion", component: Conclusion },
		{ name: "Exercise", component: Exercise01 }
	],
	files: [
		"chebi_15377_xrefs_UniProt.csv",
		"chebi_17245_xrefs_UniProt.csv",
		"chebi_27732_xrefs_UniProt.csv",
		"chebi_30050_xrefs_UniProt.csv",
		"localcurl.sh"
	],
	init: `mv /usr/local/bin/curl /usr/local/bin/curl.ori; chmod u+x localcurl.sh; cp localcurl.sh /usr/local/bin/curl;`
};
