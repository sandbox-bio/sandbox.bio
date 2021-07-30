// Bedtools
import BedtoolsIntro from "./tutorials/bedtools-intro/Intro.svelte";
import BedtoolsStep1 from "./tutorials/bedtools-intro/Step1.svelte";
import BedtoolsStep2 from "./tutorials/bedtools-intro/Step2.svelte";
import BedtoolsStep3 from "./tutorials/bedtools-intro/Step3.svelte";
import BedtoolsStep4 from "./tutorials/bedtools-intro/Step4.svelte";
import BedtoolsCpgBed from "./tutorials/bedtools-intro/data/cpg.bed";
import BedtoolsExonsBed from "./tutorials/bedtools-intro/data/exons.bed";
import BedtoolsData1Bed from "./tutorials/bedtools-intro/data/fHeart-DS15839.bed";
import BedtoolsData2Bed from "./tutorials/bedtools-intro/data/fHeart-DS16621.bed";
import BedtoolsData3Bed from "./tutorials/bedtools-intro/data/fSkin-DS19745.bed";
import BedtoolsGwasBed from "./tutorials/bedtools-intro/data/gwas.bed";
import BedtoolsChromHMMBed from "./tutorials/bedtools-intro/data/hesc.chromHmm.bed";
import BedtoolsGenomeTxt from "./tutorials/bedtools-intro/data/genome.txt";

export const config = {
	"tutorials": [
		{
			// Metadata
			"id": "bedtools-intro",
			"name": "Introduction to bedtools",
			"description": "Explore, analyze, and manipulate genomic interval <code>.bed</code> files.",
			"tools": ["bedtools"],
			"difficulty": ["beginner"],
			"author": {
				"name": "Aaron Quinlan",
				"link": "https://bioscience.utah.edu/faculty/quinlan/"
			},
			// Tutorial
			"steps": [
				{ "name": "Introduction to bedtools", "component": BedtoolsIntro },
				{ "name": "The data", "component": BedtoolsStep1 },
				{ "name": "The bedtools help", "component": BedtoolsStep2 },
				{ "name": "Bedtools 'intersect'", "component": BedtoolsStep3 },
				{ "name": "Bedtools 'intersect'", "component": BedtoolsStep4, "subtitle": "Reporting the original feature in each file" }
			],
			// Files needed
			"files": [
				{ "name": "cpg.bed", "contents": BedtoolsCpgBed },
				{ "name": "exons.bed", "contents": BedtoolsExonsBed },
				{ "name": "fHeart-DS15839.bed", "contents": BedtoolsData1Bed },
				{ "name": "fHeart-DS16621.bed", "contents": BedtoolsData2Bed },
				{ "name": "fSkin-DS19745.bed", "contents": BedtoolsData3Bed },
				{ "name": "gwas.bed", "contents": BedtoolsGwasBed },
				{ "name": "hesc.chromHmm.bed", "contents": BedtoolsChromHMMBed },
				{ "name": "genome.txt", "contents": BedtoolsGenomeTxt }
			]
		},
		{
			"id": "bowtie2-intro",
			"name": "Introduction to bowtie2",
			"description": "Align DNA sequencing reads from <code>.fastq</code> files to the Lambda phage reference genome."
		}
	]
}
