<script>
import Execute from "components/Execute.svelte";
import Exercise from "components/Exercise.svelte";

const hostname = window.location.origin;
const curl = `curl ${hostname}/data/debugging-puzzles/sequences.fa`;
const criteria = [
	{
		name: "File <code>sequences.fa</code> is not empty",
		checks: [{
			type: "file",
			path: "sequences.fa",
			action: "contents",
			commandExpected: curl
		}]
	},
	{
		name: "File <code>sequences.txt</code> contains the lines from <code>sequences.fa</code> that contain the symbol <code>&gt;</code>",
		checks: [{
			type: "file",
			path: "sequences.txt",
			action: "contents",
			commandExpected: `${curl} | grep ">"`
		}]
	}
];

const hints = [
	"The terminal interprets <code>grep > sequences.fa</code> to mean: run <code>grep</code> and output the result to <code>sequences.fa</code>, hence the file gets overwritten.",
	"Is there a way you can enclose the <code>></code> symbol so it's clear to the terminal that it is an argument to the <code>grep</code> command?"
];
</script>

Here's a fasta file that contains various sequences from the human genome and you want to know how many sequences are in the file:

<Execute command={"head sequences.fa"} />

Each FASTA record starts with a `>`, and you can assume those carets are not seen anywhere else in the record.

Naturally, you reach for good old `grep` and find lines that contain `>`. Easy peasy--click below to run the command:

<Execute command={"grep > sequences.fa"} />

That's odd, you were expecting to see lines containing `>`, not the `grep` command usage output.

Oups, we lost the data in `sequences.fa`:

<Execute command={"cat sequences.fa"} />

The file is empty, but don't worry, we can regenerate the data easily (not always the case in real life!):

<Execute command={`${curl} > sequences.fa`} />

**Your Goal**: Find the correct `grep` command that outputs lines containing `>`. Save the output as `sequences.txt`

<Exercise {criteria} {hints} />
