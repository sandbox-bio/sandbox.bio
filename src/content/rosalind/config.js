import Intro from "./steps/Intro.md";
import Exercise from "./steps/Exercise.svelte";

let rosalind = [
	{
		id: "DNA",
		title: "Counting DNA Nucleotides",
		given: "A DNA string _s_ of length at most 1000 nt.",
		return: "Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in _s_.",
		sample_data: "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
		sample_output: "20 12 17 21",
		params: ["s"]
	},
	{
		id: "RNA",
		title: "Transcribing DNA into RNA",
		given: "A DNA string _t_ having length at most 1000 nt.",
		return: "The transcribed RNA string of _t_.",
		sample_data: "GATGGAACTTGACTACGTAAATT",
		sample_output: "GAUGGAACUUGACUACGUAAAUU",
		params: ["t"]
	},
	{
		id: "REVC",
		title: "Complementing a Strand of DNA",
		given: "A DNA string _s_ of length at most 1000 bp.",
		return: "The reverse complement _s^c_ of _s_.",
		sample_data: "AAAACCCGGT",
		sample_output: "ACCGGGTTTT",
		params: ["s"]
	},
	{
		id: "FIB",
		title: "Rabbits and Recurrence Relations",
		given: "Positive integers _n <= 40_ and _k <= 5_.",
		return:
			"The total number of rabbit pairs that will be present after _n_ months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of _k_ rabbit pairs (instead of only 1 pair).",
		sample_data: "5 3",
		sample_output: "19",
		params: ["n", "k"]
	},
	{
		id: "GC",
		title: "Computing GC Content",
		given: "At most 10 DNA strings in FASTA format (of length at most 1 kbp each).",
		return:
			"The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated.",
		sample_data:
			">Rosalind_6404\nCCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC\nTCCCACTAATAATTCTGAGG\n>Rosalind_5959\nCCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT\nATATCCATTTGTCAGCAGACACGC\n>Rosalind_0808\nCCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC\nTGGGAACCTGCGGGCAGTAGGTGGAAT",
		sample_output: "Rosalind_0808\n60.919540",
		params: ["fasta"]
	},
	{
		id: "HAMM",
		title: "Counting Point Mutations",
		given: "The Hamming distance _dH(s, t)_.",
		return: "The Hamming distance _dH(s, t)_.",
		sample_data: "GAGCCTACTAACGGGAT\nCATCGTAATGACGGCCT",
		sample_output: "7",
		params: ["s", "t"]
	},
	{
		id: "IPRB",
		title: "Mendel's First Law",
		given:
			"Three positive integers _k_, _m_, and _n_, representing a population containing _k+m+n_ organisms: _k_ individuals are homozygous dominant for a factor, _m_ are heterozygous, and _n_ are homozygous recessive.",
		return:
			"The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.",
		sample_data: "2 2 2",
		sample_output: "0.78333",
		params: ["k", "m", "n"]
	},
	{
		id: "PROT",
		title: "Translating RNA into Protein",
		given: "An RNA string _s_ corresponding to a strand of mRNA (of length at most 10 kbp).",
		return: "The protein string encoded by _s_.",
		sample_data: "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA",
		sample_output: "MAMAPRTEINSTRING",
		params: ["s"]
	},
	{
		id: "SUBS",
		title: "Finding a Motif in DNA",
		given: "Two DNA strings _s_ and _t_ (each of length at most 1 kbp).",
		return: "All locations of _t_ as a substring of _s_.",
		sample_data: "GATATATGCATATACTT\nATAT",
		sample_output: "2 4 10",
		params: ["s", "t"]
	},
	{
		id: "CONS",
		title: "Consensus and Profile",
		given: "A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.",
		return:
			"A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)",
		sample_data:
			">Rosalind_1\nATCCAGCT\n>Rosalind_2\nGGGCAACT\n>Rosalind_3\nATGGATCT\n>Rosalind_4\nAAGCAACC\n>Rosalind_5\nTTGGAACT\n>Rosalind_6\nATGCCATT\n>Rosalind_7\nATGGCACT",
		sample_output: "ATGCAACT\nA: 5 1 0 0 5 5 0 0\nC: 0 0 1 4 2 0 6 1\nG: 1 1 6 3 0 1 0 0\nT: 1 5 0 0 0 1 1 6",
		params: ["fasta"]
	},
	{
		id: "FIBD",
		title: "Mortal Fibonacci Rabbits",
		given: "Positive integers _n <= 100_ and _m <= 20_.",
		return: "The total number of pairs of rabbits that will remain after the _n_-th month if all rabbits live for _m_ months.",
		sample_data: "6 3",
		sample_output: "4",
		params: ["n", "m"]
	},
	{
		id: "GRPH",
		title: "Overlap Graphs",
		given: "A collection of DNA strings in FASTA format having total length at most 10 kbp.",
		return: "The adjacency list corresponding to _O3_.  You may return edges in any order.",
		sample_data: ">Rosalind_0498\nAAATAAA\n>Rosalind_2391\nAAATTTT\n>Rosalind_2323\nTTTTCCC\n>Rosalind_0442\nAAATCCC\n>Rosalind_5013\nGGGTGGG",
		sample_output: "Rosalind_0498 Rosalind_2391\nRosalind_0498 Rosalind_0442\nRosalind_2391 Rosalind_2323",
		params: ["fasta"]
	},
	{
		id: "IEV",
		title: "Calculating Expected Offspring",
		given:
			"Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the number of couples in a population possessing each genotype pairing for a given factor. In order, the six given integers represent the number of couples having the following genotypes:",
		return:
			"The expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple has exactly two offspring.",
		sample_data: "1 0 0 1 0 1",
		sample_output: "3.5",
		params: ["a", "b", "c", "d", "e", "f"]
	},
	{
		id: "LCSM",
		title: "Finding a Shared Motif",
		given: "A collection of _k_ (_k <= 100_) DNA strings of length at most 1 kbp each in FASTA format.",
		return: "A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)",
		sample_data: ">Rosalind_1\nGATTACA\n>Rosalind_2\nTAGACCA\n>Rosalind_3\nATACA",
		sample_output: "AC",
		params: ["fasta"]
	},
	{
		id: "LIA",
		title: "Independent Alleles",
		given:
			"Two positive integers _k_ (_k <= 7_) and _N_ (_N <= 2^k_).  In this problem, we begin with Tom, who in the 0th generation has genotype Aa Bb.  Tom has two children in the 1st generation, each of whom has two children, and so on. Each organism always mates with an organism having genotype Aa Bb.",
		return:
			"The probability that at least _N_ Aa Bb organisms will belong to the _k_-th generation of Tom's family tree (don't count the Aa Bb mates at each level).  Assume that Mendel's second law holds for the factors.",
		sample_data: "2 1",
		sample_output: "0.684",
		params: ["k", "N"]
	},
	{
		id: "MPRT",
		title: "Finding a Protein Motif",
		given: "At most 15 UniProt Protein Database access IDs.",
		return:
			"For each protein possessing the N-glycosylation motif, output its given access ID followed by a list of locations in the protein string where the motif can be found.",
		sample_data: "A2Z669\nB5ZC00\nP07204_TRBM_HUMAN\nP20840_SAG1_YEAST",
		sample_output: "B5ZC00\n85 118 142 306 395\nP07204_TRBM_HUMAN\n47 115 116 382 409\nP20840_SAG1_YEAST\n79 109 135 248 306 348 364 402 485 501 614",
		params: ["ids"]
	},
	{
		id: "MRNA",
		title: "Inferring mRNA from Protein",
		given: "A protein string of length at most 1000 aa.",
		return:
			"The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)",
		sample_data: "MA",
		sample_output: "12",
		params: ["protein"]
	},
	{
		id: "ORF",
		title: "Open Reading Frames",
		given: "A DNA string _s_ of length at most 1 kbp in FASTA format.",
		return: "Every distinct candidate protein string that can be translated from ORFs of _s_. Strings can be returned in any order.",
		sample_data: ">Rosalind_99\nAGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG",
		sample_output: "MLLGSFRLIPKETLIQVAGSSPCNLS\nM\nMGMTPRLGLESLLE\nMTPRLGLESLLE",
		params: ["s"]
	},
	{
		id: "PERM",
		title: "Enumerating Gene Orders",
		given: "A positive integer _n <= 7_.",
		return: "The total number of permutations of length _n_, followed by a list of all such permutations (in any order).",
		sample_data: "3",
		sample_output: "6\n1 2 3\n1 3 2\n2 1 3\n2 3 1\n3 1 2\n3 2 1",
		params: ["n"]
	}
];

export const config = {
	id: "rosalind",
	ide: true,
	listed: false,
	name: "Rosalind Exercises",
	tags: [],
	tools: [],
	difficulty: [],
	steps: [
		{
			name: "Rosalind Exercises",
			component: Intro,
			subtitle: "Introduction",
			rosalind: {
				sample_data: "Hello!",
				sample_output: "",
				id: "welcome",
				given: "",
				return: "",
				params: ["message"]
			}
		},
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
