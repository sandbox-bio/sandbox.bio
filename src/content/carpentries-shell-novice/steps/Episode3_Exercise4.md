<script>
import Quiz from "$components/Quiz.svelte";
</script>

When run in the `alkanes` directory, which `ls` command(s) will
produce this output?

`ethane.pdb   methane.pdb`

<Quiz id="q3.4" choices={[
{ valid: false, value: "ls *t*ane.pdb"},
{ valid: false, value: "ls *t?ne.*"},
{ valid: true, value: "ls *t??ne.pdb"},
{ valid: false, value: "ls ethane.*"},
]} />
