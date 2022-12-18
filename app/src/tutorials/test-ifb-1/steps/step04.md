<script>
import Alert from "components/Alert.svelte";
import Quiz from "components/Quiz.svelte";
import Execute from "components/Execute.svelte";
</script>

The Bash shell provides you with a working space including files and directories.

A very useful command is `ls`, that **l**i**s**ts the content of a directory.
In your Unix terminal on the right, type `ls` and then press <kbd>Enter</kbd>.

The Bash shell should display `Data` that is a directory named `Data`.

Now, type the following command in your terminal (and press <kbd>Enter</kbd>) :

<Execute command="ls Data" />

The Bash shell should display the 8 files included in the `Data` directory.

Remarks:

* Pay attention to the space character between `ls` and `Data`
* Don't forget to press <kbd>Enter</kbd> to run commands
* `ls` is the command **name**
* `Data` is a directory name and an **argument** of the `ls` command

### Options

Options modify the way in which a command works.
In Bash, shell options start with a simple or double dash (`-` or `--`).

For example, we can display the size of the files using the `--size` option of the `ls` command. 
Lets try and type: 

<Execute command="ls --size Data" />

Now, the 8 files are displayed with their respective sizes (in blocks). 

You can use a short form for this option by replacing `--size` (long form) by `-s` (short form).

Usually we also use the `-h` option to display sizes in more **h**uman readable formats (_e.g._ 1K, 234M, 2G). 

You can use several options in the same command.

You can merge short form options using a single dash as prefix.
Example: type the following command in your terminal:

<Execute command="ls -sh Data" />

<Quiz id="q1" choices={[
	{ valid: true, value: "ls -s -h Data"},
	{ valid: true, value: "ls -sh Data"},
	{ valid: false, value: "ls -size -h Data"},
	{ valid: true, value: "ls --size -h Data"},
	{ valid: false, value: "ls --sizeh Data"},
	{ valid: false, value: "ls --size-h Data"},
	{ valid: true, value: "ls -h -s Data"},
	{ valid: true, value: "ls -hs Data"},
	{ valid: false, value: "ls -hsize Data"},
]}>
	<span slot="prompt">
		Among the following commands, which ones are correct?
	</span>
</Quiz>

<Alert>
### Summary

Here we learn that:
- We can execute a command by typing its name (options and/or arguments) and pressing <kbd>Enter</kbd>.
- We can add options to modify the command behavior.
- Options start with a dash `-` (short form) or a double dash `--` (long form).
- Several options of the short form can be combined (without space and after a single dash).
- Command name can be followed (or not) by an argument.
</Alert>
