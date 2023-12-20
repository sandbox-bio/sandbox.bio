<script>
import Quiz from "$components/Quiz.svelte";
import Execute from "$components/Execute.svelte";
</script>

Go to the folder `exercise-data/alkanes/`:

<Execute command="cd ~/tutorial/exercise-data/alkanes/" />

In the `alkanes` directory, imagine you have a shell script called `script.sh` containing the following commands:

```
head -n $2 $1
tail -n $3 $1
```

While you are in the `alkanes` directory, you type the following command:

```
bash script.sh '*.pdb' 1 1
```

<Quiz id="q7.1" choices={[
{ valid: false, value: `All of the lines between the first and the last lines of each file ending in .pdb in the alkanes directory`},
{ valid: true, value: `The first and the last line of each file ending in .pdb in the alkanes directory`},
{ valid: false, value: `The first and the last line of each file in the alkanes directory`},
{ valid: false, value: `An error because of the quotes around *.pdb`},
]}>
<span slot="prompt">
Which of the following outputs would you expect to see?
</span>
</Quiz>
