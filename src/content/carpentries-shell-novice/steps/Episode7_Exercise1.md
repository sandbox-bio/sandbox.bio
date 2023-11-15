<script>
import Quiz from "$components/Quiz.svelte";
import Execute from "$components/Execute.svelte";
</script>

Go to the folder `exercise-data/writing/`:

<Execute command="cd ~/tutorial/exercise-data/writing/" />

<Quiz id="q7.1" choices={[
{ valid: false, value: `grep "of" haiku.txt`},
{ valid: false, value: `grep -E "of" haiku.txt`},
{ valid: true, value: `grep -w "of" haiku.txt`},
{ valid: false, value: `grep -i "of" haiku.txt`},
]}>
<span slot="prompt">
Which command would result in the following output:
```
and the presence of absence:
```
</span>
</Quiz>
