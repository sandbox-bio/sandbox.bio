import Intro from "./steps/Intro.md";
import Exercise from "./steps/Exercise.svelte";

const rosalind = [
	{"id": "DNA", "title": "Counting DNA Nucleotides", "given": "A DNA string _s_ of length at most 1000 nt.", "return": "Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in _s_.", "sample_data": "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC", "sample_output": "20 12 17 21"}, 
	{"id": "RNA", "title": "Transcribing DNA into RNA", "given": "A DNA string _t_ having length at most 1000 nt.", "return": "The transcribed RNA string of _t_.", "sample_data": "GATGGAACTTGACTACGTAAATT", "sample_output": "GAUGGAACUUGACUACGUAAAUU"}, 
	{"id": "REVC", "title": "Complementing a Strand of DNA", "given": "A DNA string _s_ of length at most 1000 bp.", "return": "The reverse complement _s^c_ of _s_.", "sample_data": "AAAACCCGGT", "sample_output": "ACCGGGTTTT"}, 
	{"id": "FIB", "title": "Rabbits and Recurrence Relations", "given": "Positive integers _n <= 40_ and _k <= 5_.", "return": "The total number of rabbit pairs that will be present after _n_ months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of _k_ rabbit pairs (instead of only 1 pair).", "sample_data": "5 3", "sample_output": "19"}, 
	{"id": "GC", "title": "Computing GC Content", "given": "At most 10 DNA strings in FASTA format (of length at most 1 kbp each).", "return": "The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.", "sample_data": ">Rosalind_6404\nCCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC\nTCCCACTAATAATTCTGAGG\n>Rosalind_5959\nCCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT\nATATCCATTTGTCAGCAGACACGC\n>Rosalind_0808\nCCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC\nTGGGAACCTGCGGGCAGTAGGTGGAAT", "sample_output": "Rosalind_0808\n60.919540"}, 
	{"id": "HAMM", "title": "Counting Point Mutations", "given": "Two DNA strings _s_ and _t_ of equal length (not exceeding 1 kbp).", "return": "The Hamming distance _d_H(s, t)_.", "sample_data": "GAGCCTACTAACGGGAT\nCATCGTAATGACGGCCT", "sample_output": "7"}, 
	{"id": "IPRB", "title": "Mendel's First Law", "given": "Three positive integers _k_, _m_, and _n_, representing a population containing _k+m+n_ organisms: _k_ individuals are homozygous dominant for a factor, _m_ are heterozygous, and _n_ are homozygous recessive.", "return": "The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.", "sample_data": "2 2 2", "sample_output": "0.78333"}, 
	{"id": "PROT", "title": "Translating RNA into Protein", "given": "An RNA string _s_ corresponding to a strand of mRNA (of length at most 10 kbp).", "return": "The protein string encoded by _s_.", "sample_data": "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA", "sample_output": "MAMAPRTEINSTRING"}, 
	{"id": "SUBS", "title": "Finding a Motif in DNA", "given": "Two DNA strings _s_ and _t_ (each of length at most 1 kbp).", "return": "All locations of _t_ as a substring of _s_.", "sample_data": "GATATATGCATATACTT\nATAT", "sample_output": "2 4 10"}, 
	{"id": "CONS", "title": "Consensus and Profile", "given": "A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.", "return": "A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)", "sample_data": ">Rosalind_1\nATCCAGCT\n>Rosalind_2\nGGGCAACT\n>Rosalind_3\nATGGATCT\n>Rosalind_4\nAAGCAACC\n>Rosalind_5\nTTGGAACT\n>Rosalind_6\nATGCCATT\n>Rosalind_7\nATGGCACT", "sample_output": "ATGCAACT\nA: 5 1 0 0 5 5 0 0\nC: 0 0 1 4 2 0 6 1\nG: 1 1 6 3 0 1 0 0\nT: 1 5 0 0 0 1 1 6"}];

export const config = {
	id: "rosalind",
	terminal: false,
	ide: true,
	listed: false,
	name: "Rosalind Exercises",
	tags: [],
	tools: [],
	difficulty: [],
	steps: [
		{ name: "Rosalind Exercises", component: Intro, subtitle: "Introduction" },
		...rosalind.map((exercise, i) => {
			return {
				name: "Rosalind Exercises",
				component: Exercise,
				subtitle: exercise.title,
				header: i == 0,
				rosalind: exercise
			};
		})
	]
};
