import { readable } from "svelte/store";
import { config as terminalIntro } from "tutorials/terminal-basics/config.js";
import { config as igvIntro } from "tutorials/igv-intro/config.js";
import { config as bedtoolsIntro } from "tutorials/bedtools-intro/config.js";
import { config as bowtie2Intro } from "tutorials/bowtie2-intro/config.js";
import { config as samtoolsIntro } from "tutorials/samtools-intro/config.js";
import { config as jqIntro } from "tutorials/jq-intro/config.js";
import { config as awkIntro } from "tutorials/awk-intro/config.js";
import { config as fastpIntro } from "tutorials/fastp-intro/config.js";
import { config as dnaSecrets } from "tutorials/dna-secrets/config.js";
import { config as viralAmplicon } from "tutorials/viral-amplicon/config.js";
import { config as debuggingPuzzles } from "tutorials/debugging-puzzles/config.js";
import { config as playground } from "tutorials/playground/config.js";
import { config as rosalind } from "tutorials/rosalind/config.js";
import { config as IFB1 } from "tutorials/ifb-1/config.js";
import { config as IFB2 } from "tutorials/ifb-2/config.js";
// AUTO-IMPORT-HERE

// All tutorials
export const tutorials = readable([
	// Playground
	playground,
	rosalind,
	// Terminal Tutorials
	terminalIntro,
	jqIntro,
	awkIntro,
	// Bioinformatics
	bedtoolsIntro,
	bowtie2Intro,
	samtoolsIntro,
	fastpIntro,
	igvIntro,
	dnaSecrets,
	viralAmplicon,
	debuggingPuzzles,
	// IFB Tutorials (unlisted)
	IFB1,
	IFB2,
	// Add tutorials here
]);

// Playgrounds
export const playgrounds = readable([
	{
		name: "Command Line",
		description: "Command line for open-ended exploration",
		url: "/playground",
		tags: ["command line", "terminal", "bash"]
	},
	{
		name: "Awk",
		description: "Filter and wrangle tabular data",
		url: "/playgrounds?id=awk",
		tags: ["awk"]
	},
	{
		name: "Jq",
		description: "Filter and wrangle JSON data",
		url: "/playgrounds?id=jq",
		tags: ["jq", "json"]
	},
	{
		name: "Grep",
		description: "Search and filter utility",
		url: "/playgrounds?id=grep",
		tags: ["grep"]
	},
	{
		name: "Sed",
		description: "Search and replace utility",
		url: "/playgrounds?id=sed",
		tags: ["sed"]
	}
]);

// Linkouts
export const explore = readable([
	{
		name: "Bioinformatics Algorithms",
		description: "Try your hand at solving Rosalind exercises using Python.",
		url: "/rosalind",
		tags: ["rosalind", "python", "exercises"]
	},
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
