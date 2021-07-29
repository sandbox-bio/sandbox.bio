import Intro from "./tutorials/bedtools-intro/Intro.svelte";
import Step1 from "./tutorials/bedtools-intro/Step1.svelte";

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
				{ "component": Intro },
				{ "component": Step1 }
			]
		},
		{
			"id": "bowtie2-intro",
			"name": "Introduction to bowtie2",
			"description": "Align DNA sequencing reads from <code>.fastq</code> files to the Lambda phage reference genome."
		}
	]
}
