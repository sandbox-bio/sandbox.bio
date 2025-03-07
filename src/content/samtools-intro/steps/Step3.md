<script>
import Execute from "$components/Execute.svelte";
</script>

To do anything meaningful with alignment data from BWA or other aligners (which produce text-based SAM output), we need to first convert the SAM to its binary counterpart, BAM format. The binary format is much easier for computer programs to work with. However, it is consequently very difficult for humans to read. More on that later.

To convert SAM to BAM, we use the `samtools view` command. We must specify that we want the output to be BAM (by default it produces SAM) with the `-b` option. Samtools follows the UNIX convention of sending its output to the UNIX STDOUT, so we need to use `-o` to create a BAM file from the output.

<Execute command="samtools view -b sample.sam -o sample.bam" />

Now, we can use the `samtools view` command to convert the BAM to SAM so we mere mortals can read it:

<Execute command="samtools view sample.bam | head -n 5" />
