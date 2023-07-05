<script>
import Execute from "$components/Execute.svelte";
import Exercise from "$components/Exercise.svelte";

const criteria = [{
	name: "File <code>alignments.fixed.bam</code> can be queried at a random position",
	checks: [{
		type: "file",
		path: "alignments.fixed.bam",
		action: "contents",
		commandObserved: `samtools view alignments.fixed.bam ref2:10-11`,
		commandExpected: `samtools sort alignments.bam -o /shared/tmp/__debuggingpuzzles.bam; samtools index /shared/tmp/__debuggingpuzzles.bam; samtools view /shared/tmp/__debuggingpuzzles.bam ref2:10-11`
	}]
}];

const hints = [
	"The first step is to create an index for <code>alignments.bam</code>. Check out <code>samtools --help</code> for the command you need to create an index.",
	"Indexing using <code>samtools index alignments.bam</code> fails because the alignments are not sorted, so we first need to sort the file!",
	"Look into using <code>samtools sort</code> to sort <code>alignments.bam</code> before indexing the resulting file."
];
</script>

You just finished aligning sequencing data using your favorite read aligner, which output the file `alignments.bam`:

<Execute command="samtools view alignments.bam" />

Although you can view the entire file, you want to extract alignments to specific genomic coordinates.

But you get an error when you run the following `samtools` command:

<Execute command="samtools view alignments.bam ref2:10-11" />

**Your Goal**: Create a file `alignments.fixed.bam` from the alignments in `alignments.bam` so that the following query does not give an error:

<Execute command="samtools view alignments.fixed.bam ref2:10-11" />

<Exercise {criteria} {hints} />
