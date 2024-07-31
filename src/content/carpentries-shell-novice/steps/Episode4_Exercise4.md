<script>
import Quiz from "$components/Quiz.svelte";
</script>

Consider the file `exercise-data/animal-counts/animals.csv`:

```
2012-11-05,deer,5
2012-11-05,rabbit,22
2012-11-05,raccoon,7
2012-11-06,rabbit,19
...
```

The `uniq` command has a `-c` option which gives a count of the number of times a line occurs in its input.

<Quiz id="q4.4" choices={[
{ valid: false, value: "sort animals.csv | uniq -c"},
{ valid: false, value: "sort -t, -k2,2 animals.csv | uniq -c"},
{ valid: false, value: "cut -d, -f 2 animals.csv | uniq -c"},
{ valid: true, value: "cut -d, -f 2 animals.csv | sort | uniq -c"},
{ valid: false, value: "cut -d, -f 2 animals.csv | sort | uniq -c | wc -l"},
]}>
<span slot="prompt">
    What command would you use to produce a table that shows the total count of each type of animal in the file?
</span>
</Quiz>
