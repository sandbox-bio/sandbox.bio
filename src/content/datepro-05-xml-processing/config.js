// Steps
import Intro from "./steps/Intro.md";
import Step01 from "./steps/Step01.md";
import Step02 from "./steps/Step02.md";
import Step03 from "./steps/Step03.md";
import Conclusion from "./steps/Conclusion.md";
import Exercise01 from "./exercises/Exercise01.md";

export const config = {
	id: "datepro-05-xml-processing",
	pwd: "datepro-05-xml-processing",
	name: "Section 3.7: XML Processing",
	subtitle: `by <a href="https://webpages.ciencias.ulisboa.pt/~fjcouto/" target="_blank">Francisco M. Couto</a>`,
	description: "Searching, Filtering, and Extracting XML Data",
	tags: ["proteins", "search", "XPath query", "XML", "UniProt"],
	tools: ["xmllint", "grep", "expr", "sort", "wc", "cut"],
	difficulty: ["intermediate"],
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "Human proteins filter", component: Step01 },
		{ name: "PubMed identifiers extraction", component: Step02 },
		{ name: "Complex elements", component: Step03 },
		{ name: "Conclusion", component: Conclusion },
		{ name: "Exercise", component: Exercise01 }
	],
	files: [
		"chebi_27732_A2AGL3.xml",
		"chebi_27732_B0LPN4.xml",
		"chebi_27732_E9PZQ0.xml",
		"chebi_27732_E9Q401.xml",
		"chebi_27732_F1LMY4.xml",
		"chebi_27732_O14174.xml",
		"chebi_27732_O94236.xml",
		"chebi_27732_O94513.xml",
		"chebi_27732_P11716.xml",
		"chebi_27732_P21817.xml",
		"chebi_27732_P38689.xml",
		"chebi_27732_P54867.xml",
		"chebi_27732_Q15413.xml",
		"chebi_27732_Q5AHG6.xml",
		"chebi_27732_Q8N490.xml",
		"chebi_27732_Q92736.xml",
		"chebi_27732_Q9TS33.xml",
		"chebi_27732_Q9VD76.xml",
		"chebi_27732_Q9VKA5.xml",
		"chebi_27732_Q9VSH2.xml"
	]
};
