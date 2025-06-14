<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

Now, let's begin the phylogenetic analysis by performing Multiple Sequence Alignment (MSA) for our SARS-CoV-2 sequences. Recall that like human genomes, viral genomes can "evolve" (i.e. mutate) as they replicate. Common mutations include substitutions, insertions, and deletions. Given a set of viral sequences, each of which differ from the sequence of the common ancestor by a series of mutations, it is our job to first "line up" each of our sequences so that each position in the "alignment" of our sequences corresponds to the same position in the sequence of the common ancestor. 

We will use [ViralMSA](https://github.com/niemasd/ViralMSA) to perform MSA. Check out how ViralMSA should be used with <Execute command="ViralMSA.py -h" inline />.

To run MSA for the sequences in `sarscov2_sequences.fas`:

<Execute command={"ViralMSA.py \\ --omit_ref \\ -s sarscov2_sequences.fas \\ -r sarscov2_reference.fas \\ -o ViralMSA_Out"} />

Take a quick look at `./ViralMSA_Out/sarscov2_sequences.fas.aln`. This file is still in FASTA format, but what changed?

<Quiz
	id="step1-quiz3"
	choices={[
		{ valid: false, value: `There are more sequences in the file` },
		{ valid: false, value: `There are fewer sequences in the file` },
		{ valid: false, value: `The sequence identifiers have been truncated` },
		{ valid: true, value: `Dashes are inserted to represent gaps in the alignment` },
    ]}>
	<span slot="prompt"></span>
</Quiz>
