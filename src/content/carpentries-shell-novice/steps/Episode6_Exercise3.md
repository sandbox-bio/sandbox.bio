<script>
import Exercise from "$components/Exercise.svelte";
import Execute from "$components/Execute.svelte";

const criteria = [
{
	name: "Script <code>longest.sh</code> exists",
	checks: [{
		type: "file",
		path: "exercise-data/longest.sh",
		action: "exists"
	}]
},
{
	name: "The output file <code>longest.txt</code> exists: <code>bash longest.sh alkanes pdb > longest.txt</code>",
	checks: [{
		type: "file",
		path: "exercise-data/longest.txt",
		action: "exists"
	}]
},
{
	name: "The file <code>longest.txt</code> contains each file's unique species",
	checks: [{
		type: "file",
		path: "exercise-data/longest.txt",
		action: "contents",
		commandExpected: `cd exercise-data/ && wc -l alkanes/*.pdb | sort -n | tail -n 2 | head -n 1`
	}]
}];
</script>

Go to the folder `exercise-data/`:

<Execute command="cd $TUTORIAL/exercise-data/" />

Write a shell script called `longest.sh` that takes the name of a directory and a filename extension as its arguments, and prints out the name of the file with the most lines in that directory with that extension. For example:

```bash
bash longest.sh alkanes pdb
```

would print the name of the `.pdb` file in `alkanes/` that has the most lines.

Feel free to test your script on another directory e.g.

```bash
bash longest.sh writing txt
```

<Exercise {criteria} />
