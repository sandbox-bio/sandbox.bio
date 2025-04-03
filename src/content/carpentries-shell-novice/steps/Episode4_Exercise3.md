<script>
import Alert from "$components/Alert.svelte";
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

For the file `animals.csv`, consider the following command:

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
