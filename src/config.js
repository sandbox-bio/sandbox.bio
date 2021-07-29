// import bla from "./1-intro-to-bedtools/data/cpg.bed"
// import bedtools_1_text from "./1-intro-to-bedtools/step-1.md"

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
					"name": "Synopsis",
					"text": "yes"
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
