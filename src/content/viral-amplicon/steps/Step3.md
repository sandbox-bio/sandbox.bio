<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

Now that we have a SAM file, we'll use `viral_consensus` to perform consensus sequence calling using our mapped reads. Internally, `viral_consensus` first trims the reads by performing two types of trimming:

- **Primer Trimming:** In order to sequence specific target regions of interest, the amplicon sequencing protocol results in the addition of a <Link href="https://en.wikipedia.org/wiki/Primer_(molecular_biology)#Uses_of_synthetic_primers">synthetic primer</Link> at the beginning/end of a read. In "Primer Trimming," we remove these primers from the beginning/end of our reads
- **Quality Trimming:** Like all technologies, sequencing machines are error-prone, and each nucleotide contained within a read has a <Link href="https://en.wikipedia.org/wiki/FASTQ_format#Quality">quality score</Link> associated with it. In "Quality Trimming," we remove low-quality bases from the beginning/end of our reads

Then, `viral_consensus` uses the trimmed reads to count the nucleotides at each position in order to call a "consensus sequence": a viral genome sequence containing the "consensus" (i.e., most abundant) nucleotide at each position. Both steps (trimming and consensus sequence calling) happen in a single `viral_consensus` command. Run the following:

<Execute command={"viral_consensus -i reads.mapped.sam \ -r $REF_FASTA \ -o consensus.fa \ -op position_counts.tsv \ -oi insertion_counts.json \ -p $PRIMER_BED \ -po 5"} />

Let's break this seemingly complex command into its individual components to make some sense of it:

- `-i reads.mapped.sam` tells `viral_consensus` that our input file is `reads.mapped.sam`
- `-r $REF_FASTA` tells `viral_consensus` that we want to use the file `$REF_FASTA` as our reference genome FASTA
- `-o consensus.fa` tells `viral_consensus` that we want to write the output consensus genome sequence to the file `consensus.fa`
  - The output file is a FASTA file
- `-op position_counts.tsv` tells `viral_consensus` to write the position-by-position nucleotide counts to the file `position_counts.tsv`
  - This is optional, but the nucleotide counts can be useful for quality-checking the resulting consensus sequence
- `-oi insertion_counts.json` tells `viral_consensus` to write the insertion counts to the file `insertion_counts.json`
  - This is optional, but the insertion counts can be useful for quality-checking the resulting consensus sequence
- `-p $PRIMER_BED` tells `viral_consensus` that our primer BED file is `$PRIMER_BED`
  - A <Link href="https://en.wikipedia.org/wiki/BED_(file_format)">BED file</Link> contains a list of genomic regions
  - The primer BED file contains the the start and end positions (with respect to the reference genome) of our amplicon sequencing primers
- `-po 5` tells `viral_consensus` to trim any reads within 5 positions of a primer
  - In other words, instead of strictly trimming only reads that start within a primer, reads that fall outside of a primer but are close enough (5 positions) will be included in Primer Trimming

After running the above command, we will have successfully called a consensus genome sequence and written the resulting sequence to the file `consensus.fa`.

To see the contents of the consensus genome sequence output file, run the following:

<Execute command="cat consensus.fa" />

This consensus genome FASTA file is the other key result of our amplicon sequencing analysis, and it is what will be used in any downstream analyses (e.g. transmission clustering, phylogenetic inference, lineage assignment, etc.).

We will have also output the nucleotide counts (`position_counts.tsv`) and insertion counts (`insertion_counts.json`), which we can use for quality-checking. Specifically, we can see how many reads supported each nucleotide or insertion in the consensus sequence that we called.

To see the contents of the position counts TSV file, run the following:

<Execute command="head -5 position_counts.tsv" />

To see the contents of the insertion counts JSON file, run the following:

<Execute command="cat insertion_counts.json" />
