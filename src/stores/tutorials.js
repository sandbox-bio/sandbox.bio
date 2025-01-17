import { readable } from "svelte/store";
import { config as _template } from "$content/_template/config.js";
import { config as terminalIntro } from "$content/terminal-basics/config.js";
import { config as igvIntro } from "$content/igv-intro/config.js";
import { config as bedtoolsIntro } from "$content/bedtools-intro/config.js";
import { config as bowtie2Intro } from "$content/bowtie2-intro/config.js";
import { config as samtoolsIntro } from "$content/samtools-intro/config.js";
import { config as jqIntro } from "$content/jq-intro/config.js";
import { config as awkIntro } from "$content/awk-intro/config.js";
import { config as fastpIntro } from "$content/fastp-intro/config.js";
import { config as dnaSecrets } from "$content/dna-secrets/config.js";
import { config as viralAmplicon } from "$content/viral-amplicon/config.js";
import { config as debuggingPuzzles } from "$content/debugging-puzzles/config.js";
import { config as playground } from "$content/playground/config.js";
import { config as ifblinuxbasics1 } from "$content/ifb-linux-basics-1/config.js";
import { config as ifblinuxbasics2 } from "$content/ifb-linux-basics-2/config.js";
import { config as ifblinuxbasics3 } from "$content/ifb-linux-basics-3/config.js";
import { config as carpentriesShellNovice } from "$content/carpentries-shell-novice/config";
import { config as blastIntro } from "$content/blast-intro/config";
import { config as jellyfishIntro } from "$content/jellyfish-intro/config";
import { config as seqkitIntro } from "$content/seqkit-intro/config";
import { config as datepro1 } from "$content/datepro-01-unix-shell/config";
import { config as datepro2 } from "$content/datepro-02-data-retrieval/config";
import { config as datepro3 } from "$content/datepro-03-data-extraction/config";
import { config as datepro4 } from "$content/datepro-04-task-repetition/config";
import { config as datepro5 } from "$content/datepro-05-xml-processing/config";
import { config as datepro6 } from "$content/datepro-06-text-retrieval/config";
import { config as datepro7 } from "$content/datepro-07-text-processing/config";
import { config as datepro8 } from "$content/datepro-08-semantic-processing/config";
import { config as ska2Intro } from "$content/ska2-intro/config";
import { env } from "$env/dynamic/public";

// All tutorials
export const tutorials = readable([
	// Playground
	playground,
	// Terminal Tutorials
	terminalIntro,
	carpentriesShellNovice,
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
	jellyfishIntro,
	blastIntro,
	seqkitIntro,
	ska2Intro,
	// Community tutorials
	ifblinuxbasics1,
	ifblinuxbasics2,
	ifblinuxbasics3,
	datepro1,
	datepro2,
	datepro3,
	datepro4,
	datepro5,
	datepro6,
	datepro7,
	datepro8,
	// Template tutorial
	_template
]);

// Tutorial listings
export const categories = readable([
	{
		name: "Recently added",
		icon: "star-fill",
		tutorials: (env?.PUBLIC_USE_PRD_ASSETS ? [_template] : []).concat([ska2Intro, seqkitIntro, carpentriesShellNovice, blastIntro]),
		mailinglist: true
	},
	{
		name: "Data exploration",
		icon: "compass-fill",
		tutorials: [terminalIntro, carpentriesShellNovice, jqIntro, awkIntro]
	},
	{
		name: "File formats",
		icon: "file-earmark-binary-fill",
		tutorials: [bedtoolsIntro, samtoolsIntro, seqkitIntro]
	},
	{
		name: "Quality control",
		icon: "bookmark-check-fill",
		tutorials: [fastpIntro, igvIntro]
	},
	{
		name: "Data analysis",
		icon: "cpu-fill",
		tutorials: [bowtie2Intro, blastIntro, jellyfishIntro, dnaSecrets, viralAmplicon, ska2Intro, debuggingPuzzles]
	}
]);

// Playgrounds
export const playgrounds = readable([
	{
		name: "Command Line",
		description: "Command line for open-ended exploration",
		url: "/tutorials/playground",
		tags: ["command line", "terminal", "bash"]
	},
	{
		name: "Awk",
		description: "Filter and wrangle tabular data",
		url: "/playgrounds/awk",
		tags: ["awk"]
	},
	{
		name: "Jq",
		description: "Filter and wrangle JSON data",
		url: "/playgrounds/jq",
		tags: ["jq", "json"]
	},
	{
		name: "Grep",
		description: "Search and filter utility",
		url: "/playgrounds/grep",
		tags: ["grep", "regex"]
	},
	{
		name: "Sed",
		description: "Search and replace utility",
		url: "/playgrounds/sed",
		tags: ["sed", "regex"]
	}
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
