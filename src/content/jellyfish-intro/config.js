// Steps
import Intro from "./steps/Intro.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";
import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "jellyfish-intro",
	name: "K-mer counting with Jellyfish",
	icon: "list-ol",
	subtitle: `by <a href="https://robert.bio" target="_blank">Robert Aboukhalil</a>`,
	description: "Learn the basics of k-mer counting using Jellyfish",
	tags: ["kmer", "counting", "jellyfish"],
	tools: ["ls", "jellyfish", "seqtk", "wgsim"],
	difficulty: ["beginner"],
	new: true,
	steps: [
		{ name: "K-mer counting with Jellyfish", component: Intro },
		{ name: "Count k-mers in genomes", component: Step1 },
		{ name: "Exercise: K-mer counting", component: Step2 },
		{ name: "Query k-mers", component: Step3 },
		{ name: "Count k-mers in sequencing reads", component: Step4 },
		{ name: "The end", component: Conclusion }
	],
	files: ["dengue.fa", "chikungunya.fa"]
};
