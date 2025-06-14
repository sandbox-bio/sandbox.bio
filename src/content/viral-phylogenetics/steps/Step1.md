<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
import Quiz from "$components/Quiz.svelte";
</script>

For this tutorial, we will be using SARS-CoV-2 whole genome sequences collected from samples from real people! These sequences can be found in the file `sarscov2_sequences.fas`. Let's start by taking a peek at this file.

Try <Execute command="less -S sarscov2_sequences.fas" inline /> to take a look at the SARS-CoV-2 sequences we will be using. To exit the view, you can simply type `q` for `quit`.

`sarscov2_sequences.fas` is in FASTA format, the standardized file format for nucleotide sequence data.

What information do the characters after the `>` represent? Take a look at the FASTA format <Link href="https://www.ncbi.nlm.nih.gov/genbank/fastaformat/">documentation</Link>! 

<Quiz
	id="step1-quiz1"
	choices={[
		{ valid: false, value: `the length of the sequence` },
		{ valid: true, value: `the unique sequence identifier` },
		{ valid: false, value: `the sequence quality` },
		{ valid: false, value: `the random sequence barcode` },
    ]}>
	<span slot="prompt"></span>
</Quiz>

Next, let's figure out how many sequences are in `sarscov2_sequences.fas`. To do so, we can use the `grep "<pattern>" <file>` command, which will enable us to search for a particular text pattern in a file. For example, `grep "abc" test.txt` will return all lines containing the string `"abc"`. Use the `grep` command to determine how many sequences are in `sarscov2_sequences.fas`.

Hint: add `"^"` to the beginning of `"pattern"` to limit the search to lines that begin with `"pattern"`. The `wc -l` command may also be helpful!

How many sequences are in `sarscov2_sequences.fas`?

<Quiz
	id="step1-quiz2"
	choices={[
		{ valid: false, value: `1` },
		{ valid: false, value: `9` },
		{ valid: false, value: `100` },
		{ valid: true, value: `11` },
    ]}>
	<span slot="prompt"></span>
</Quiz>

If you're up for an extra challenge, try to determine the length of each sequence in `sarscov2_sequences.fas`!
