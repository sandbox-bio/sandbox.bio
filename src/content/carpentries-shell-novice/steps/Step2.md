<script>
import Quiz from "$components/Quiz.svelte";
import Execute from "$components/Execute.svelte";
import Exercise from "$components/Exercise.svelte";
</script>

Suppose we have several hundred genome data files named `basilisk.dat`, `minotaur.dat`, and
`unicorn.dat`.
For this example, we'll use the `exercise-data/creatures` directory which only has three
example files,
but the principles can be applied to many many more files at once.

<Execute command="cd ~/tutorial/exercise-data/creatures" />

The structure of these files is the same: the common name, classification, and updated date are
presented on the first three lines, with DNA sequences on the following lines.

Let's look at the files:

<Execute command="head -n 5 basilisk.dat minotaur.dat unicorn.dat" />

We would like to print out the classification for each species, which is given on the second
line of each file.
For each file, we would need to execute the command `head -n 2` and pipe this to `tail -n 1`.
We'll use a loop to solve this problem, but first let's look at the general form of a loop,
using the pseudo-code below:

```bash
# The word "for" indicates the start of a
# "For-loop" command
for thing in list_of_things
# The word "do" indicates the start of job
# execution list
do
  # Indentation within the loop is not required,
  # but aids legibility
  operation_using/command $thing
# The word "done" indicates the end of a loop
done
```

and we can apply this to our example like this:

<Execute command="for filename in basilisk.dat minotaur.dat unicorn.dat \ndo\n  echo $filename; \n  head -n 2 $filename | tail -n 1; \ndone" />
