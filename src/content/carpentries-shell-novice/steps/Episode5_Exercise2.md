<script>
import Quiz from "$components/Quiz.svelte";
</script>

What would be the output of running the following loop in the `exercise-data/alkanes` directory?

```bash
for filename in c*
do
    ls $filename
done
```

<Quiz id="q5.2" choices={[
{ valid: false, value: "No files are listed."},
{ valid: false, value: "All files are listed."},
{ valid: false, value: "Only cubane.pdb, octane.pdb and pentane.pdb are listed."},
{ valid: true, value: "Only cubane.pdb is listed."},
]}>
<span slot="prompt">
</span>
</Quiz>

How would the output differ from using this command instead (`*c*` instead of `c*`):

```bash
for filename in *c*
do
    ls $filename
done
```

<Quiz id="q5.2b" choices={[
{ valid: false, value: "The same files would be listed."},
{ valid: false, value: "All the files are listed this time."},
{ valid: false, value: "No files are listed this time."},
{ valid: true, value: "The files cubane.pdb and octane.pdb will be listed."},
{ valid: false, value: "Only the file octane.pdb will be listed."},
]}>
<span slot="prompt">
</span>
</Quiz>
