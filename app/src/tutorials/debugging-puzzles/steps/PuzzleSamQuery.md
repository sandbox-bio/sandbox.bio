<script>
import Execute from "components/Execute.svelte";
import Exercise from "components/Exercise.svelte";

const criteria = [{
	name: "Converted <code>alignments.bam</code>",
	checks: [{
		type: "file",
		path: "exons.bed",
		action: "contents",
		commandExpected: `sed 's/ /\t/g' exons.bed`
	}]
}];

const hints = [];
</script>

You have just finished aligning sequencing data using your favorite read aligner, which output the file `alignments.bam`:

<Execute command="samtools view alignments.bam" />

Although you can view the entire file, you want to extract alignments to specific genomic coordinates.

But you get an error when you run the following `samtools` command:

<Execute command="samtools view alignments.bam ref2:10-11" />

**Your Goal**: Create a file `fixed.bam` from the alignments in `alignment.bam` so that the following query does not give an error:

<Execute command="samtools view fixed.bam ref2:10-11" />

<Exercise {criteria} {hints} />
