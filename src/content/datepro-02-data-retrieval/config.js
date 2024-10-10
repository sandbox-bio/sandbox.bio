// Steps
import Intro from "./steps/Intro.md";
import Step01 from "./steps/Step01.md";
import Step02 from "./steps/Step02.md";
import Step03 from "./steps/Step03.md";
import Conclusion from "./steps/Conclusion.md";
import Exercise01 from "./exercises/Exercise01.md";

export const config = {
	id: "datepro-01-data-retrieval",
	pwd: "datepro-01-data-retrieval",
	name: "Data and Text Processing: Data Retrieval",
	subtitle: `by <a href="https://webpages.ciencias.ulisboa.pt/~fjcouto/" target="_blank">Francisco M. Couto</a>`,
	description: "Data Retrieval",
	tags: ["curl", "redirection", "parameters", "CSV"],
	tools: ["curl", "cat", "nano", "chmod"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Introduction", component: Intro },
		{ name: "Web Identifiers", component: Step01 },
		{ name: "Data Retrieval", component: Step02 },
		{ name: "Manage Output", component: Step03 },
		{ name: "Conclusion", component: Conclusion },
		{ name: "Exercise", component: Exercise01 },
	],
	files: [
		"chebi_15377_xrefs_UniProt.csv",
        "chebi_17245_xrefs_UniProt.csv",
        "chebi_27732_xrefs_UniProt.csv",
        "chebi_27732_xrefs_UniProt.xls",
        "chebi_27732_xrefs_UniProt.xml",
        "chebi_30050_xrefs_UniProt.csv",
        "localcurl.sh",
		],
		init: `alias curl='data/localcurl.sh'; chmod u+x data/localcurl.sh`
};
