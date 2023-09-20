<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

Next, we'll use `ivar` to trim the reads in the sorted BAM file we created in the previous step. Specifically, `ivar` will perform two types of trimming:

- **Primer Trimming:** In order to sequence specific target regions of interest, the amplicon sequencing protocol results in the addition of a <Link href="https://en.wikipedia.org/wiki/Primer_(molecular_biology)#Uses_of_synthetic_primers">synthetic primer</Link> at the beginning/end of a read. In "Primer Trimming," we remove these primers from the beginning/end of our reads
- **Quality Trimming:** Like all technologies, sequencing machines are error-prone, and each nucleotide contained within a read has a <Link href="https://en.wikipedia.org/wiki/FASTQ_format#Quality">quality score</Link> associated with it. In "Quality Trimming," we remove low-quality bases from the beginning/end of our reads

Run the following:

<Execute command="ivar trim -x 5 -e \ -i untrimmed.sorted.bam \ -b $PRIMER_BED \ -p trimmed.unsorted" />

Let's break this seemingly complex command into its individual components to make some sense of it:

- `trim` tells `ivar` that we want to use its trimming functionality
- `-x 5` tells `ivar` to trim any reads within 5 positions of a primer
  - In other words, instead of strictly trimming only reads that start within a primer, reads that fall outside of a primer but are close enough (5 positions) will be included in Primer Trimming
- `-e` tells `ivar` to include reads that don't have primers
  - By default, `ivar trim` excludes reads that don't have primers
  - In other words, `-e` tells `ivar` to trim _all_ reads
- `-i untrimmed.sorted.bam` tells `ivar` that our input file is `untrimmed.sorted.bam`
- `-b $PRIMER_BED` tells `ivar` that our primer BED file is `$PRIMER_BED`
  - A <Link href="https://en.wikipedia.org/wiki/BED_(file_format)">BED file</Link> contains a list of genomic regions
  - The primer BED file contains the the start and end positions (with respect to the reference genome) of our amplicon sequencing primers
- `-p trimmed.unsorted` tells `ivar` that we want the prefix of our output file (i.e., the file name without the file extension) to be `trimmed.sorted`
  - `ivar` will automatically add the `.bam` file extension to the prefix we specify
  - In other words, this will result in an output file named `trimmed.unsorted.bam`

After running the above command, we will have successfully trimmed our mapped reads and written the trimmed (no longer sorted) results to the file `trimmed.unsorted.bam`.

To see the first few lines of the (unsorted) trimmed BAM output file in the human-readable SAM format, run the following:

<Execute command="samtools view -h trimmed.unsorted.bam | \ head -n 5" />
