// Steps
import Intro from "./steps/Intro.md";
import Step01 from "./steps/Step01.md";
import Step02 from "./steps/Step02.md";
import Step03 from "./steps/Step03.md";
import Conclusion from "./steps/Conclusion.md";
import Exercise01 from "./exercises/Exercise01.md";

export const config = {
	id: "datepro-02-data-retrieval",
	next: "datepro-03-data-extraction",
	pwd: "datepro-02-data-retrieval",
	name: "Sections 3.3 and 3.4: Data Retrieval",
	subtitle: `by <a href="https://webpages.ciencias.ulisboa.pt/~fjcouto/" target="_blank">Francisco M. Couto</a>`,
	description: "Utilizing web services to access the ChEBI database, retrieving proteins associated with specific compounds",
	tags: ["curl", "redirection", "parameters", "CSV"],
	tools: ["curl", "cat", "nano", "chmod"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "Web Identifiers", component: Step01 },
		{ name: "Data Retrieval", component: Step02 },
		{ name: "Manage Output", component: Step03 },
		{ name: "Exercise", component: Exercise01 },
		{ name: "Conclusion", component: Conclusion },
	],
	files: [
		"chebi_15377_xrefs_UniProt.csv",
		"chebi_17245_xrefs_UniProt.csv",
		"chebi_27732_xrefs_UniProt.csv",
		"chebi_27732_xrefs_UniProt.xls",
		"chebi_27732_xrefs_UniProt.xml",
		"chebi_30050_xrefs_UniProt.csv",
		"localcurl.sh"
	],
	// init: "chmod +x localcurl.sh && alias curl='./localcurl.sh'"
	init: `mv /usr/local/bin/curl /usr/local/bin/curl.ori; chmod u+x localcurl.sh; cp localcurl.sh /usr/local/bin/curl;`
};
