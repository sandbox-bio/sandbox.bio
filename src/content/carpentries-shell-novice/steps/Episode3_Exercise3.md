<script>
import Quiz from "$components/Quiz.svelte";
</script>

What is the output of the closing `ls` command in the sequence shown below?

```bash
$ pwd
/Users/jamie/data

$ ls
proteins.dat

$ mkdir recombined
$ mv proteins.dat recombined/
$ cp recombined/proteins.dat ../proteins-saved.dat
$ ls
```

<Quiz id="q3.3" choices={[
{ valid: false, value: "proteins-saved.dat recombined"},
{ valid: true, value: "recombined"},
{ valid: false, value: "proteins.dat recombined"},
{ valid: false, value: "proteins-saved.dat"},
]}>
<span slot="prompt">
What is the output of the `ls` command above?
</span>
</Quiz>
