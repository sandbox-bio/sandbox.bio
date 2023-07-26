<script>
import Quiz from "components/Quiz.svelte";
import Execute from "components/Execute.svelte";
</script>

Using a Bash shell, you will be able to write and execute Unix commands. 
A Unix command is made of distinct parts: a **command name** and, if needed, **options** and **arguments**.

The space character is mandatory to separate all these elements (command name, options and arguments). 
You may now understand why a space in a file or a directory name is a very bad idea within an Unix environment.

Be very careful while writing your Unix commands: lowercases and uppercases are not the same! 
Unix commands are usually in lowercase.

Now, type the following command in your terminal and then press <kbd>Enter</kbd> key:

<Execute command="date" />

<Quiz id="q1" choices={[
	{ valid: true, value: "yes"},
	{ valid: false, value: "no"},
]}>
	<span slot="prompt">
		Does the terminal display the current date?
	</span>
</Quiz>
