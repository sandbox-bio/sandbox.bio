import { readable } from "svelte/store";
import { config as terminalIntro } from "tutorials/terminal-basics/config.js";
import { config as bedtoolsIntro } from "tutorials/bedtools-intro/config.js";
import { config as bowtie2Intro } from "tutorials/bowtie2-intro/config.js";
import { config as samtoolsIntro } from "tutorials/samtools-intro/config.js";
import { config as jqIntro } from "tutorials/jq-intro/config.js";
import { config as dnaSecrets } from "tutorials/dna-secrets/config.js";
import { config as playground } from "tutorials/playground/config.js";
import { config as rosalind } from "tutorials/rosalind/config.js";

// All tutorials
export const tutorials = readable([
	// Playground
	playground,
	rosalind,
	// Terminal Tutorials
	terminalIntro,
	jqIntro,
	// Bioinformatics
	bedtoolsIntro,
	bowtie2Intro,
	samtoolsIntro,
	dnaSecrets
]);

// Linkouts
export const explore = readable([
	{
		name: "Align DNA sequences",
		description: "Explore the Smith-Waterman and Needleman-Wunsch sequence alignment algorithms.",
		url: "https://alignment.sandbox.bio/",
		tags: ["smith-waterman", "needleman-wunsch"]
	},
	{
		name: "Simulate DNA sequences",
		description: "Simulate DNA sequencing reads with <code>wgsim</code>.",
		url: "https://wgsim.sandbox.bio/",
		tags: ["wgsim"]
	},
	{
		name: "t-SNE algorithm",
		description: "Run t-SNE on single-cell sequencing data.",
		url: "https://tsne.sandbox.bio/",
		tags: ["t-SNE"]
	},
	{
		name: "QC reports for FASTQ files",
		description: "Generate data quality reports with <code>fastp</code>.",
		url: "https://fastq.sandbox.bio/",
		tags: ["fastp"]
	}
]);