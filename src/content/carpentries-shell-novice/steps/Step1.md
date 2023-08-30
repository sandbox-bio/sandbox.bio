<script>
import Quiz from "$components/Quiz.svelte";
import Execute from "$components/Execute.svelte";
import Exercise from "$components/Exercise.svelte";

let criteria = [
{
	name: "File <code>unique.csv</code> exists",
	checks: [{
		type: "file",
		path: "exercise-data/animal-counts/unique.csv",
		action: "exists"
	}]
},
{
	name: "File <code>unique.csv</code> contains the unique rows from the file <code>animals.csv</code>",
	checks: [{
		type: "file",
		path: "exercise-data/animal-counts/unique.csv",
		action: "contents",
		commandExpected: `cut -d , -f 2 exercise-data/animal-counts/animals.csv | sort | uniq`
	}]
}];
</script>

Now that we know a few basic commands,
we can finally look at the shell's most powerful feature:
the ease with which it lets us combine existing programs in new ways.
We'll start with the directory `exercise-data/alkanes`
that contains six files describing some simple organic molecules.

<Execute command="cd exercise-data/alkanes" />

The `.pdb` extension indicates that these files are in Protein Data Bank format,
a simple text format that specifies the type and position of each atom in the molecule.

<Execute command="ls" />

Let's run an example command:

<Execute command="wc cubane.pdb" />

`wc` is the 'word count' command:
it counts the number of lines, words, and characters in files (returning the values
in that order from left to right).

If we run the command `wc *.pdb`, the `*` in `*.pdb` matches zero or more characters,
so the shell turns `*.pdb` into a list of all `.pdb` files in the current directory:

<Execute command="wc *.pdb" />

Note that `wc *.pdb` also shows the total number of all lines in the last line of the output.

If we run `wc -l` instead of just `wc`,
the output shows only the number of lines per file:

<Execute command="wc -l *.pdb" />

The `-m` and `-w` options can also be used with the `wc` command to show
only the number of characters or the number of words, respectively.

[...]

## Combining multiple commands

Nothing prevents us from chaining pipes consecutively.
We can for example send the output of `wc` directly to `sort`,
and then send the resulting output to `head`.
This removes the need for any intermediate files.

We'll start by using a pipe to send the output of `wc` to `sort`:

<Execute command="wc -l *.pdb | sort -n" />

We can then send that output through another pipe, to `head`, so that the full pipeline becomes:

<Execute command="wc -l *.pdb | sort -n | head -n 1" />

This is exactly like a mathematician nesting functions like _log(3x)_
and saying 'the log of three times _x_'.
In our case,
the algorithm is 'head of sort of line count of `*.pdb`'.

[...]

## Piping Commands Together

In our current directory, we want to find the 3 files which have the least number of
lines.

<Quiz id="q1" choices={[
{ valid: false, value: "wc -l * > sort -n > head -n 3"},
{ valid: false, value: "wc -l * | sort -n | head -n 1-3"},
{ valid: false, value: "wc -l * | head -n 3 | sort -n"},
{ valid: true, value: "wc -l * | sort -n | head -n 3"},
]}>
<span slot="prompt">
Which command listed below would work?
</span>
</Quiz>

[...]

## Pipe Construction

For the file `animals.csv` from the previous exercise, consider the following command:

<Execute command="cut -d , -f 2 ~/tutorial/exercise-data/animal-counts/animals.csv" />

The `cut` command is used to remove or 'cut out' certain sections of each line in the file,
and `cut` expects the lines to be separated into columns by a <kbd>Tab</kbd> character.
A character used in this way is a called a **delimiter**.
In the example above we use the `-d` option to specify the comma as our delimiter character.
We have also used the `-f` option to specify that we want to extract the second field (column).

The `uniq` command filters out adjacent matching lines in a file.
How could you extend this pipeline (using `uniq` and another command) to find
out what animals the file contains (without any duplicates in their
names)?

Store the result in unique.csv:

<Exercise {criteria} />
