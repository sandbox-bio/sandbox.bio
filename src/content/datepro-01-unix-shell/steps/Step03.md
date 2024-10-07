<script>
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
</script>

#### Line breaks
A Unix file represents a single line break by a line feed character, instead of
two characters (carriage return and line feed) used by Windows. So, if we
are using a text editor in Windows, we must be careful to use one that lets
us save it as Unix file, for example the open source Notepad++. If we are
using a text editor in macOS we also need to be careful in saving it in text
format. In case we do not have such text editor, we can also remove the extra
carriage return by using the command line tool `tr`, that replaces and deletes
characters:

<Execute command="tr -d '\r' < reversemyfile.sh > reversemyfilenew.sh" />

The `-d` option of `tr` is used to remove a given character from the input, in
this case tr will delete all carriage returns (`r`). Many command line options
can be used in short form using a single dash (`-`), or in a long form using two
dashes (`--`). In this tool, using the `--delete` option is equivalent to the `-d`
option. Long forms are more self-explanatory, but they take longer to type
and occupy more space. 

#### Redirection operator

The `>` character represents a redirection operator that moves the results
being displayed at the standard output (our terminal) to a given file. The `<`
character represents a redirection operator that works on the opposite direction, i.e. opens a given file and uses it as the standard input.
We should note that `cat` received the filename as an input argument,
while `tr` can only receive the contents of the file through the standard input.

Instead of providing the filename as argument, the cat command can also
receive the contents of a file through the standard input, and produce the
same output:

<Execute command="cat < myfile.txt" />

The previous `tr` command used a new file for the standard output, because we cannot use the same file to read and write at the same time. To keep the same filename, we have to move the new file by using the mv command:

<Execute command="mv reversemyfilenew.sh reversemyfile.sh" />

And restore the execution permissions:

<Execute command="chmod u+x reversemyfile.sh" />

#### Save output

We can now save the output into another file named mynewfile.txt by typing:

<Execute command="./reversemyfile.sh myfile.txt > mynewfile.txt" />

Again, to check if the file was really created, we can use the `cat` tool:

<Execute command="cat mynewfile.txt" />

Or, we can reverse it again by typing:

<Execute command="./reversemyfile.sh mynewfile.txt" />

Of course, the result should exactly be the original contents of `myfile.txt`.


#### Debug

If something is not working well, we can debug the entire script by typing:

<Execute command="bash -x reversemyfile.sh myfile.txt" />

Our terminal will not only display the resulting text, but also the command
line tools executed preceded by the plus character (`+`)
Alternatively, we can add the `set -x` command line in our script to start the
debugging mode, and `set +x` to stop it.
