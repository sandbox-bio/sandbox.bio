<script>
import Exercise from "$components/Exercise.svelte";
import Execute from "$components/Execute.svelte";

const hints = [
    "Use <code>man grep</code> to look for how to grep text recursively in a directory",
    "Use <code>man cut</code> to select more than one field in a line"
];
const criteria = [
{
	name: "Script <code>count-species.sh</code> exists",
	checks: [{
		type: "file",
		path: "exercise-data/count-species.sh",
		action: "exists"
	}]
},
{
	name: "Calling the script creates the file <code>bear.txt</code>: <code>bash count-species.sh bear animal-counts/</code> ",
	checks: [{
		type: "file",
		path: "exercise-data/bear.txt",
		action: "exists"
	}]
},
{
	name: "The file <code>bear.txt</code> contains the list of dates and the number of bears seen",
	checks: [{
		type: "file",
		path: "exercise-data/bear.txt",
		action: "contents",
		commandExpected: `grep -w bear -r exercise-data/animal-counts/ | cut -d : -f 2 | cut -d , -f 1,3`
	}]
}];
</script>

Go to the folder `exercise-data/`:

<Execute command="cd ~/tutorial/exercise-data/" />

Leah has several hundred data files saved in one directory, each of which is formatted like this:

```js
2012-11-05,deer,5
2012-11-05,rabbit,22
2012-11-05,raccoon,7
2012-11-06,rabbit,19
2012-11-06,deer,2
2012-11-06,fox,4
2012-11-07,rabbit,16
2012-11-07,bear,1
```

She wants to **write a shell script that takes a species as the first command-line argument and a directory as the second argument**. The script should return one file called `<species>.txt` containing a list of dates and the number of that species seen on each date. For example using the data shown above, `rabbit.txt` would contain:

```js
2012-11-05,22
2012-11-06,19
2012-11-07,16
```

Below, each line contains an individual command, or pipe. Arrange their sequence in one command in order to achieve Leah’s goal:

```bash
cut -d : -f 2
>
|
grep -w $1 -r $2
|
$1.txt
cut -d , -f 1,3
```

<Exercise {criteria} {hints} />
