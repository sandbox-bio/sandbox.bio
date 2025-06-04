<script>
import Exercise from "$components/Exercise.svelte";
import Execute from "$components/Execute.svelte";

const criteria = [
{
	name: "Script <code>species.sh</code> exists",
	checks: [{
		type: "file",
		path: "exercise-data/species.sh",
		action: "exists"
	}]
},
{
	name: "The output file <code>unique.txt</code> exists: <code>bash species.sh animal-counts/animals.csv > unique.txt</code>",
	checks: [{
		type: "file",
		path: "exercise-data/unique.txt",
		action: "exists"
	}]
},
{
	name: "The file <code>unique.txt</code> contains each file's unique species",
	checks: [{
		type: "file",
		path: "exercise-data/unique.txt",
		action: "contents",
		commandExpected: `cut -d , -f 2 exercise-data/animal-counts/animals.csv | sort | uniq`
	}]
}];
</script>

Go to the folder `exercise-data/`:

<Execute command="cd $TUTORIAL/exercise-data/" />

Leah has several hundred data files, each of which is formatted like this:

```js
2013-11-05,deer,5
2013-11-05,rabbit,22
2013-11-05,raccoon,7
2013-11-06,rabbit,19
2013-11-06,deer,2
2013-11-06,fox,1
2013-11-07,rabbit,18
2013-11-07,bear,1
```

An example of this type of file is given in `animal-counts/animals.csv`.

We can use the command `cut -d , -f 2 animal-counts/animals.csv | sort | uniq` to produce the unique species in `animals.csv`. In order to avoid having to type out this series of commands every time, a scientist may choose to write a shell script instead.

Write a shell script called `species.sh` that takes any number of filenames as command-line arguments and uses a variation of the above command to print a list of the unique species appearing in each of those files separately.

<Exercise {criteria} />
