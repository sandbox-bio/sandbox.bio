<script>
import Alert from "$components/Alert.svelte";
import Quiz from "$components/Quiz.svelte";
import Execute from "$components/Execute.svelte";
</script>

A shell is a software program that interprets and executes command lines
given by the user in consecutive lines of text. A shell script is a list of such
command lines. The command line usually starts by invoking a command
line tool.

Unix shell was developed to manage Unix-like operating systems, but due to their usefulness nowadays they are available is most personal computers using Linux, macOS or Windows operating systems.

First we need to know what is our current directory.
Type: <Execute command={"pwd"} inline />.

After hitting enter, the command will show the full path of the directory
(folder) of our computer in which the shell is working on. The dollar sign in
the left is only to indicate that this is a command to be executed directly in
the shell. A curved arrow can appear in the right each time a command does
not fit in the available width of a page, and has to be presented in multiple
lines


To understand a command line tool, such as pwd, we can type man followed by the name of the tool.

For example, we can type <Execute command={"man pwd"} inline />. to learn
more about pwd (do not forget to hit enter, and press q to quit).

We can also
learn more about man by typing man man. A shorter alternative to man, is to
add the --help option after any command tool. For example, we can type
pwd --help to have a more concise description of pwd.

Type <Execute command={"ls"} inline /> in the command line.

