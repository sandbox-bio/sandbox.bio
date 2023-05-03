// Steps
import Intro from "./steps/Intro.md";
import Quiz from "./steps/Quiz.md";
// import Step2 from "./steps/Step2.md";
// import Step3 from "./steps/Step3.md";
// import Step4 from "./steps/Step4.md";
// import Step5 from "./steps/Step5.md";
// import Step6 from "./steps/Step6.md";
// import Step7 from "./steps/Step7.md";
// import Step8 from "./steps/Step8.md";
// import Conclusion from "./steps/Conclusion.md";

export const config = {
	id: "igv-intro",
	pwd: "igv-intro",
	name: "Visualize variants with IGV",
	subtitle: `by <a href="https://rnabio.org" target="_blank">The Griffith Lab</a>`,
	description: "Distinguish variants from artifacts.",
	tags: ["igv"],
	tools: [],
	difficulty: ["beginner"],
	igv: true,
	igvConfig: {
		locus: "all",
		tracks: [{
			url: "https://rnabio.org/assets/module_2/HCC1143.normal.21.19M-20M.bam",
			indexURL: "https://rnabio.org/assets/module_2/HCC1143.normal.21.19M-20M.bam.bai",
		}],
	},
	steps: [
		{ name: "Visualize variants with IGV", component: Intro },
		{ name: "Quiz", component: Quiz },
		// { name: "Filter data", component: Step1, subtitle: "Select Elements", header: true },
		// { name: "Filter data", component: Step2, subtitle: "Select Arrays" },
		// { name: "Filter data", component: Step3, subtitle: "Putting Elements in an Array" },
		// { name: "Filter data", component: Step4, subtitle: "Select Multiple Fields" },
		// { name: "Filter data", component: Step5, subtitle: "Putting Elements Into an Object" },
		// { name: "Summarize data", component: Step6, subtitle: "Sorting and Counting", header: true },
		// { name: "Summarize data", component: Step7, subtitle: "Pipes and Filters" },
		// { name: "Summarize data", component: Step8, subtitle: "Maps and Selects" },
		// { name: "The end", component: Conclusion, subtitle: "In Review", header: true }
	],
	files: [],
};

/*

Get familiar with the interface

## Navigation

You should see listing of chromosomes in this reference genome. Choose 1, for chromosome 1.

### Chromosome chooser

Navigate to chr1:10,000-11,000 by entering this into the location field (in the top-left corner of the interface) and clicking Go.

This shows a window of chromosome 1 that is 1,000 base pairs wide and beginning at position 10,000.

### Colors

IGV displays the sequence of letters in a genome as a sequence of colours (e.g. A = green, C = blue, etc.).

This makes repetitive sequences, like the ones found at the start of this region, easy to identify.

Zoom in a bit more using the + button to see the individual bases of the reference genome sequence.

### Gene model

Genes are represented as lines and boxes. Lines represent intronic regions, and boxes represent exonic regions.

The arrows indicate the direction/strand of transcription for the gene. When an exon box become narrower in height, this indicates a UTR.

When loaded, tracks are stacked on top of each other. You can identify which track is which by consulting the label to the left of each track.


## Loading Read Alignments

We will be using the breast cancer cell line HCC1143 to visualize alignments. For speed, only a small portion of chr21 will be loaded (19M:20M).

Copy the files to your local drive, and in IGV choose File > Load from File..., select the bam file, and click OK. Note that the bam and index files must be in the same directory for IGV to load these properly.




## Visualizing read alignments

Navigate to a narrow window on chromosome 21: `chr21:19,480,041-19,480,386`.

// Experiment with the various settings by right clicking the read alignment track and toggling the options.
// Think about which would be best for specific tasks (e.g. quality control, SNP calling, CNV finding).

### Colors

You will see reads represented by grey or white bars stacked on top of each other, where they were aligned to the reference genome.

The reads are pointed to indicate their orientation (i.e. the strand on which they are mapped).

Click on any read and notice that a lot of information is available.

Once you select a read, you will learn what many of these metrics mean, and how to use them to assess the quality of your datasets.

At each base that the read sequence **mismatches** the reference, the colour of the base represents the letter that exists in the read (using the same colour legend used for displaying the reference).



## Inspecting SNPs, SNVs, and SVs

In this section we will be looking in detail at 8 positions in the genome, and determining whether they represent real events or artifacts.

### Two neighbouring SNPs

Navigate to region `chr21:19,479,237-19,479,814`.

Note two heterozygous variants, one corresponds to a known dbSNP (`G/T` on the right) the other does not (`C/T` on the left)

Zoom in and center on the `C/T` SNV on the left, sort by base (window `chr21:19,479,321` is the SNV position)

Sort alignments by `base`

Color alignments by `read strand`

#### Good quality SNVs/SNPs

Notes:

* High base qualities in all reads except one (where the alt allele is the last base of the read)
* Good mapping quality of reads, no strand bias, allele frequency consistent with heterozygous mutation

Question(s):
* What does *Shade base by quality* do? How might this be helpful?
* How does Color by *read strand* help?

### Homopolymer region with indel

Navigate to position `chr21:19,518,412-19,518,497`

#### Example 2a
Group alignments by read strand

Center on the A within the homopolymer run (chr21:19,518,470), and Sort alignments by -> base

// image: poor baseQ support in fwd reads
// https://rnabio.org/assets/module_2/example2a.png

#### Example 2b

Center on the one base deletion (chr21:19,518,452), and Sort alignments by -> base

// image: alternating insertions and deletions in reverse strand reads
// https://rnabio.org/assets/module_2/example2b.png

Notes:

* The alt allele is either a deletion or insertion of one or two Ts
* The remaining bases are mismatched, because the alignment is now out of sync
* The dpSNP entry at this location (rs74604068) is an A->T, and in all likelihood an artifact
* i.e. the common variants from dbSNP include some cases that are actually common misalignments caused by repeats




### Coverage by GC

Navigate to position chr21:19,611,925-19,631,555. Note that the range contains areas where coverage drops to zero in a few places.

* Use Collapsed view
* Use Color alignments by -> insert size and pair orientation
* Load GC track
* See concordance of coverage with GC content

Question:

* Why are there blue and red reads throughout the alignments?


### Heterozygous SNPs on different alleles

Navigate to region chr21:19,666,833-19,667,007

Sort by base (at position chr21:19,666,901)

Note:
* There is no linkage between alleles for these two SNPs because reads covering both only contain one or the other


### Low mapping quality

Navigate to region chr21:19,800,320-19,818,162

Load repeat track (File -> Load from server...)

Notes:
* Mapping quality plunges in all reads (white instead of grey). Once we load repeat elements, we see that there are two LINE elements that cause this.



### Homozygous deletion

Navigate to region chr21:19,324,469-19,331,468

* Turn on View as Pairs and Expanded view
* Use Color alignments by -> insert size and pair orientation
* Sort reads by insert size
* Click on a red read pair to pull up information on alignments

Notes:
* Typical insert size of read pair in the vicinity: 350bp
* Insert size of red read pairs: 2,875bp
* This corresponds to a homozygous deletion of 2.5kb


### Mis-alignment

Navigate to region chr21:19,102,154-19,103,108

Notes:
* This is a position where AluY element causes mis-alignment.
* Misaligned reads have mismatches to the reference and well-aligned reads have partners on other chromosomes where additional ALuY elements are encoded.
* Zoom out until you can clearly see the contrast between the difficult alignment region (corresponding to an AluY) and regions with clean alignments on either side


### Translocation

Navigate to region chr21:19,089,694-19,095,362

* Expanded view
* Group alignments by -> pair orientation
* Color alignments by -> pair orientation


Notes:
* Many reads with mismatches to reference
* Read pairs in RL pattern (instead of LR pattern)
* Region is flanked by reads with poor mapping quality (white instead of grey)
* Presence of reads with pairs on other chromosomes (coloured reads at the bottom when scrolling down)




*/