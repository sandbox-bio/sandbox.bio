<script>
import Quiz from "$components/Quiz.svelte";
</script>

Also in the `exercise-data/alkanes` directory, what would be the output of the following loop?

```bash
for datafile in *.pdb
do
    cat $datafile >> all.pdb
done
```

<Quiz id="q5.4" choices={[
{ valid: false, value: "All of the text from cubane.pdb, ethane.pdb, methane.pdb, octane.pdb, and pentane.pdb would be concatenated and saved to a file called all.pdb."},
{ valid: false, value: "The text from ethane.pdb will be saved to a file called all.pdb."},
{ valid: true, value: "All of the text from cubane.pdb, ethane.pdb, methane.pdb, octane.pdb, pentane.pdb and propane.pdb would be concatenated and saved to a file called all.pdb."},
{ valid: false, value: "All of the text from cubane.pdb, ethane.pdb, methane.pdb, octane.pdb, pentane.pdb and propane.pdb would be printed to the screen and saved to a file called all.pdb."},
]} />
