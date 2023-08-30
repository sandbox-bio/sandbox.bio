<script>
import Quiz from "$components/Quiz.svelte";
import Execute from "$components/Execute.svelte";
import Exercise from "$components/Exercise.svelte";
</script>

[...]

Let's start by going back to `alkanes/`:
<Execute command="cd ~/tutorial/exercise-data/alkanes" />

and creating a new file, `middle.sh` which will
become our shell script:

<Execute command="nano middle.sh" />

The command `nano middle.sh` opens the file `middle.sh` within the text editor 'nano'
(which runs within the shell).
If the file does not exist, it will be created.
We can use the text editor to directly edit the file by inserting the following line:

<Execute command="head -n 15 octane.pdb | tail -n 5" />

This is a variation on the pipe we constructed earlier, which selects lines 11-15 of
the file `octane.pdb`. Remember, we are _not_ running it as a command just yet;
we are only incorporating the commands in a file.

Then we save the file (`Ctrl-O` in nano) and exit the text editor (`Ctrl-X` in nano).
Check that the directory `alkanes` now contains a file called `middle.sh`.

Once we have saved the file,
we can ask the shell to execute the commands it contains.
Our shell is called `bash`, so we run the following command:

<Execute command="bash middle.sh" />

Sure enough,
our script's output is exactly what we would get if we ran that pipeline directly.
