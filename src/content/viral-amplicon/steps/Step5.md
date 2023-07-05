<script>
import Execute from "$components/Execute.svelte";
</script>

Now that we have a BAM file that is trimmed and sorted, we'll use `samtools` to generate a pile-up file that lists all of the nucleotides and corresponding quality scores covering each position of the reference genome. Run the following:

<Execute command="samtools mpileup \ -A -aa -d 0 -Q 0 \ --reference $REF_FASTA \ trimmed.sorted.bam > pileup.txt" />

Let's break this seemingly complex command into its individual components to make some sense of it:

- First, we're calling the following command: `samtools mpileup -A -aa -d 0 -Q 0 --reference $REF_FASTA trimmed.sorted.bam`
  - `mpileup` tells `samtools` that we want to use its **m**ultiway **pileup** functionality
  - `-A` tells `samtools` not to discard "**a**nomalous" read pairs (i.e., "orphans"), which are reads that mapped but whose pair failed to map
  - `-aa` tells `samtools` to output **a**bsolutely **a**ll positions, including reference positions with 0 depth or unused reference sequences
  - `-d 0` tells `samtools` to not use a "max **d**epth" limit
    - In other words, even if a position has *millions* of reads spanning a given position, we want all nucleotides at that position to appear in the output
  - `-Q 0` tells `samtools` not to perform any base **q**uality filtering
    - Normally, you can use `-Q THRESHOLD` to have `samtools` ignore any nucleotides with a base quality less than `THRESHOLD`
    - We've already accounted for *some* low-quality base calls when we performed Quality Trimming previously
    - We'll handle ignoring any remaining low-quality base calls when we perform our downstream analyses
  - `--reference $REF_FASTA` tells `samtools` that we want to use the file `$REF_FASTA` as our reference genome FASTA
  - Lastly, we're specifying the sorted BAM file from which we want to compute a pile-up file: `trimmed.sorted.bam`
- Then, we're using the redirect character `>` to redirect the output of `samtools` (via standard output) to a file: `pileup.txt`

After running the above command, we will have sucessfully computed a pile-up file and written the results to the file `pileup.txt`.

To see the first few lines of the pile-up output file, run the following:

<Execute command="head -n 5 pileup.txt" />
