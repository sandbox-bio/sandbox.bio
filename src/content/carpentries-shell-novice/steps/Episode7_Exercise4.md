<script>
import Quiz from "$components/Quiz.svelte";
import Execute from "$components/Execute.svelte";
</script>

Go to the folder `exercise-data/`:

<Execute command="cd ~/tutorial/exercise-data/" />

The `-v` option to `grep` inverts pattern matching, so that only lines which do not match the pattern are printed. Given that, which of the following commands will find all `.dat` files in `creatures` except `unicorn.dat`? Once you have thought about your answer, you can test the commands in the `exercise-data/` directory.

<Quiz id="q7.4" choices={[
{ valid: true, value: `find creatures -name "*.dat" | grep -v unicorn`},
{ valid: false, value: `find creatures -name *.dat | grep -v unicorn`},
{ valid: false, value: `grep -v "unicorn" $(find creatures -name "*.dat")`},
{ valid: false, value: `None of the above.`},
]}>
<span slot="prompt">
Which of the following is correct?
</span>
</Quiz>
