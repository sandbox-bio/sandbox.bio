// Steps
import Intro from "./steps/Intro.md";
import Step0 from "./steps/Step0.md";
import Step1 from "./steps/Step1.md";
import Step2 from "./steps/Step2.md";
import Step3 from "./steps/Step3.md";
import Step4 from "./steps/Step4.md";

export const config = {
	id: "ska2-intro",
	name: "Building trees with SKA",
	subtitle: `by <a href="https://github.com/vrbouza" target="_blank">Víctor Rodríguez Bouza</a> and <a href="https://maklin.fi/" target="_blank">Tommi Mäklin</a>`,
	description: "Use ska.rust to compare and align closely related small genomes using split k-mers",
	tags: ["ska2", "k-mers"],
	tools: ["ska"],
	difficulty: ["beginner"],
	steps: [
		{ name: "Building trees with SKA", component: Intro },
		{ name: "Why SKA is good for building phylogenetic trees", component: Step0 },
		{ name: "Indexing the assemblies", component: Step1 },
		{ name: "Producing the SNP alignment", component: Step2 },
		{ name: "Building the trees", component: Step3 },
		{ name: "Working with references", component: Step4 }
	],
	files: [
		"create_tree.py",
		"assemblies/GCA_000005845.2.fna.gz",
		"assemblies/LA160.fa.gz",
		"assemblies/LA189.fa.gz",
		"assemblies/LA207.fa.gz",
		"assemblies/LA243.fa.gz",
		"assemblies/ska_input_ref.tsv",
		"assemblies/ska_input.tsv"
	]
};
