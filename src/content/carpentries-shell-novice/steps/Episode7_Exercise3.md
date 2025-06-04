<script>
import Exercise from "$components/Exercise.svelte";
import Execute from "$components/Execute.svelte";

const hints = [
	"One solution might employ the commands <code>grep</code> and <code>wc</code> and a <code>|</code>, while another might utilize <code>grep</code> options. There is often more than one way to solve a programming task, so a particular solution is usually chosen based on a combination of yielding the correct result, elegance, readability, and speed."
];
const criteria = [
{
	name: "File <code>stats.txt</code> exists",
	checks: [{
		type: "file",
		path: "exercise-data/writing/stats.txt",
		action: "exists"
	}]
},
{
	name: "The file <code>stats.txt</code> contains how often each name was mentioned",
	checks: [{
		type: "file",
		path: "exercise-data/writing/stats.txt",
		action: "contents",
		commandExpected: `for sis in Jo Meg Beth Amy; do echo "$sis:"; grep -ow $sis exercise-data/writing/LittleWomen.txt | wc -l; done`
	}]
}];
</script>

Go to the folder `exercise-data/writing/`:

<Execute command="cd $TUTORIAL/exercise-data/writing/" />

You and your friend, having just finished reading Little Women by Louisa May Alcott, are in an argument. Of the four sisters in the book, Jo, Meg, Beth, and Amy, your friend thinks that Jo was the most mentioned. You, however, are certain it was Amy. Luckily, you have a file `LittleWomen.txt` containing the full text of the novel (`exercise-data/writing/LittleWomen.txt`). Using a `for` loop, how would you tabulate the number of times each of the four sisters is mentioned?

The file `stats.txt` should match the following format:

```js
Jo:
111
Meg:
222
Beth:
333
Amy:
444
```

<Exercise {criteria} {hints} />
