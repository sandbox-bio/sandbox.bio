<script>
import Execute from "components/Execute.svelte";
</script>

Next, we'll use `samtools` to re-sort our trimmed BAM file. Run the following:

<Execute command="samtools sort -o trimmed.sorted.bam trimmed.unsorted.bam" />

This command is very similar to one that we saw previously:

- `sort` tells `samtools` that we want to use its sorting functionality
- `-o trimmed.sorted.bam` tells `samtools` that we want to write the output BAM results to the file `trimmed.sorted.bam` (rather than printing them to standard output)
  - BAM is a binary (i.e., non-human-readable) and (potentially) compressed file format that stores the same information as the SAM format
- Lastly, we're specifying the (unsorted) BAM file we want to sort: `trimmed.unsorted.bam`

After running the above command, we will have sucessfully sorted our (unsorted) trimmed mapped reads and written the sorted results to the file `trimmed.sorted.bam`.

To see the first few lines of the sorted trimmed BAM output file in the human-readable SAM format, run the following:

<Execute command="samtools view -h trimmed.sorted.bam | \ head -n 5" />
