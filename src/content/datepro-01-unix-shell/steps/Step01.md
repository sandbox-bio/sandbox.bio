<script>
import Execute from "$components/Execute.svelte";
</script>

A shell is a software program that interprets and executes command lines
given by the user in consecutive lines of text. A shell script is a list of such
command lines. The command line usually starts by invoking a command
line tool.

Unix shell was developed to manage Unix-like operating systems, but due to their usefulness nowadays they are available is most personal computers using Linux, macOS or Windows operating systems.

#### Current directory

As our first command line, we can type:

<Execute command="pwd" />

The command will show the full path of the directory (folder) in which the shell is working on.
To understand any command line tool, such as `pwd`, we can 
add the --help option after any command tool. For example, 
to have a more concise description of pwd we can type:

<Execute command="pwd --help" /> 

As our second command line, we can type:

<Execute command="ls" />

It will show the list of files in the current directory.

#### Change directory

To change the directory, we can use another command line tool, the `cd`
(change directory) followed by the new path. In a personal computer with Linux we may want to use the Documents directory. Type:

<Execute command="cd Documents" /> 

We will receive the error saying that is no such file or directory.
This happens because we do not have that directory. 
To create we can use the `mkdir` command 

<Execute command="mkdir Documents" /> 

And now try again:

<Execute command="cd Documents" /> 

Now to see what changed we can type:

<Execute command="pwd" /> 

And if we want to return to the parent directory, we only need to use the
two dots:

<Execute command="cd .." /> 

And if we want to return to the home directory, we only need to use the
tilde character (~):

<Execute command="cd ~" /> 

To return the original `pwd` output, we should type:

<Execute command="cd /root/tutorial" />

/root/tutorial

#### Useful key combinations

Every time the terminal is blocked by any reason, we can press both the
control (**Ctrl**) and **C key** at the same time . This usually cancels the current
tool being executed. For example, try using the `cd` command with only one
single quote:

<Execute command="cd '" /> 

This will block the terminal, because it is still waiting for a second single
quote that closes the argument. Now press **Ctrl-c**, and the command will be
aborted. Now we can type again the previous command, but instead of pressing
**Ctrl-c** we may also press **Ctrl-d**. The combination **Ctrl-d** indicates the terminal that it is the end of input. So, in this case, the `cd` command will not
be canceled, but instead it is executed without the second single quote and
therefore a syntax error will be shown on our display.

Other useful key combinations are the **Ctrl-l** that when pressed cleans the
terminal display,  and the **Ctrl-Insert** and **Shift-Insert** that when pressed
copy and paste the selected text, respectively.

#### Shell version
The following examples will probably work in any Unix shell, but if we want
to be certain that we are using bash we can type the following command,
and check if the output says bash:

<Execute command="ps -p $$" />

This tool shows information about active processes running in our computer. The `-p` option selects a given process, and in this case `$$` represents the process running in our terminal application. In most terminal applications bash is the default shell. If this is not our case, we may need to type bash, hit enter and now we are using bash.

Now that we know how to use a shell, the next step is to write and run a very simple script that reverses the order of the lines in a text file.
