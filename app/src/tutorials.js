import { readable } from "svelte/store";
import { config as bedtoolsIntro } from "tutorials/bedtools-intro/config.js";
import { config as bowtie2Intro } from "tutorials/bowtie2-intro/config.js";
import { config as samtoolsIntro } from "tutorials/samtools-intro/config.js";
import { config as dnaSecrets } from "tutorials/dna-secrets/config.js";

export const tutorials = readable([
	bedtoolsIntro,
	bowtie2Intro,
	samtoolsIntro,
	dnaSecrets
]);

export const explore = readable([
	{
		name: "Align DNA sequences",
		description: "Explore the Smith-Waterman and Needleman-Wunsch sequence alignment algorithms.",
		url: "https://alignment.sandbox.bio/",
		tags: ["smith-waterman", "needleman-wunsch"]
	},
	{
		name: "Simulate DNA sequences",
		description: "Use <code>wgsim</code> to simulate DNA sequencing reads.",
		url: "https://wgsim.sandbox.bio/",
		tags: ["wgsim"]
	},
	{
		name: "The t-SNE algorithm",
		description: "Run the t-SNE algorithm on single-cell sequencing data.",
		url: "https://tsne.sandbox.bio/",
		tags: ["t-SNE"]
	},
	{
		name: "QC reports for FASTQ files",
		description: "Generate data quality reports for FASTQ files using the <code>fastp</code> tool.",
		url: "https://fastq.sandbox.bio/",
		tags: ["fastp"]
	}
])
