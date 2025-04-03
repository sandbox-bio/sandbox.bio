<script>
import Execute from "$components/Execute.svelte";
</script>

#### Data file

We start by creating a file named `myfile.txt` using any text editor, for example:

<Execute command="nano myfile.txt" />

Copy and paste the following lines:

```text
line 1
line 2
line 3
line 4
```

Exit from the editor with **Ctrl-X**, and then press **Y** and **Enter** to save the file.

#### File contents

We cannot forget to save it in our working directory, and check if it has the
proper filename extension. To check if the file is really on our working directory, we can type:

<Execute command="cat myfile.txt" />

The contents of the file should appear in our terminal. `cat` is a simple
command line tool that receives a filename as argument and displays its contents on the screen.

#### Reverse file contents

An alternative to `cat` tool is the `tac` tool. To try it, we only need to type:

<Execute command="tac myfile.txt" />

The contents of the file should also appear in our terminal, but now in the
reverse order.

#### My first script

Now we can create a script file named `reversemyfile.sh` by using the text editor:

<Execute command="nano reversemyfile.sh" />

And add the following line:

```text
tac $1
```

`$1` represents the first argument after the script filename when invoking it.

Exit from the editor with **Ctrl-X**, and then press **Y** and **Enter** to save the file.

We cannot forget to save the file in our working directory. A simple way to ensure this is to check if the file and its contents are present.

<Execute command="cat reversemyfile.sh" />

We could also add the shebang `#!/bin/bash` as the first line
of the script, which would specify that it should be executed using the Bash
shell. However, for simplicity we will not use any shebang in our tutorials.

### Permissions

A script also needs permission to be executed, so every time we create a new
script file we need to type:

<Execute command="chmod u+x reversemyfile.sh" />

The command line tool `chmod` just gave the user (`u`) permissions to execute
(`+x`).

Finally, we can execute the script by providing the `myfile.txt` as argument:

<Execute command="./reversemyfile.sh myfile.txt" />

The contents of the file should appear in our terminal in the reverse order.
Congratulations, we made our first script work!

In case the terminal responds with a permission denied error message, we
need to check if the `chmod` was done correctly. In case the terminal responds with a failed to open `myfile.txt` error message, we need to check if the `reversemyfile.sh` file was saved in the same directory, as explaine before.

If we give more arguments, they will be ignored:

<Execute command="./reversemyfile.sh myfile.txt myotherfile.txt 'myother file.txt'" />

The output will be exactly the same because our script does not use `$2` and
`$3`, that in this case will represent `myotherfile.txt` and my `'other file.txt'`, respectively. We should note that when containing spaces, the argument must
be enclosed by single quotes.

The next and final step will about how to save the output of the script.
