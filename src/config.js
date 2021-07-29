// Bedtools
import BedtoolsIntro from "./tutorials/bedtools-intro/Intro.svelte";
import BedtoolsStep1 from "./tutorials/bedtools-intro/Step1.svelte";
//
import BedtoolsCpgBed from "./tutorials/bedtools-intro/data/cpg.bed";

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
				{
					"component": BedtoolsIntro
				},
				{
					"name": "The data",
					"component": BedtoolsStep1
				}
			],
			// Files needed
			"files": [
				{
					"name": "cpg.bed",
					"contents": BedtoolsCpgBed
				}
			]
		},
		{
			"id": "bowtie2-intro",
			"name": "Introduction to bowtie2",
			"description": "Align DNA sequencing reads from <code>.fastq</code> files to the Lambda phage reference genome."
		}
	]
}
