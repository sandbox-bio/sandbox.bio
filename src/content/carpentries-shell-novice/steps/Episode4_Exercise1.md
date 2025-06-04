<script>
import Quiz from "$components/Quiz.svelte";
</script>

We have already met the `head` command, which prints lines from the start of a file.
`tail` is similar, but prints lines from the end of a file instead.

Consider the file `exercise-data/animal-counts/animals.csv`.
After these commands, select the answer that
corresponds to the file `animals-subset.csv`:

```bash
$ head -n 3 animals.csv > animals-subset.csv
$ tail -n 2 animals.csv >> animals-subset.csv
```

<Quiz id="q4.1" choices={[
{ valid: false, value: "The first three lines of `animals.csv`"},
{ valid: false, value: "The last two lines of `animals.csv`"},
{ valid: true, value: "The first three lines and the last two lines of `animals.csv`"},
{ valid: false, value: "The second and third lines of `animals.csv`"},
]} />
