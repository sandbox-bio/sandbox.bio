<script>
import Alert from "components/Alert.svelte";
import Execute from "components/Execute.svelte";
</script>

<Alert>
	You can skip this tutorial if you're familiar with the command-line and are comfortable using `ls`, `head`, `tail`, `grep`, and environment variables.
</Alert>

In this tutorial, you'll get familiar with the interface on your right: **the terminal**, also known as _command-line_, _shell_, _Bash_, _UNIX command line_, etc.

Note that the sandbox.bio terminal might not behave exactly the same way it does on your computer&mdash;it is after all running inside your browser!&mdash;but you'll be able to apply what you learn here to other terminals you'll encounter in the wild.

Let's get started!

---

The terminal is a text-based interface that interprets commands and outputs the result to your screen.

The first command we'll try is `echo`, which simply returns the string you provide it:

<Execute command='echo "Hello World"' />

Try variations of that command with your own strings.

<!-- For example, add spaces between the words. Then try removing the `"` quotes, and notice the differences.

<Alert>
	A recurring theme of the command-line is that spacing and quotes are important.
</Alert> -->

Just like the graphical interface you use for exploring files on your computer, you always work within one folder when you're using the terminal.

For most of the tutorials on sandbox.bio, you won't have to worry about which folder you're in, but you can use `pwd` to **p**rint the **w**orking **d**irectory:

<Execute command="pwd" />

<!-- We can use `cd` to **c**hange the **d**irectory we're in:

<Execute command="cd /shared" />

And `pwd` will reflect that:

<Execute command="pwd" />

Let's go back to our previous directory:

<Execute command="cd /shared/data" /> -->

To list the files inside this folder, use `ls`:

<Execute command="ls" />



---

A common task we perform on the command-line is previewing files using the `head`, `tail` and `grep` commands.

To view the first 10 lines of a file, use `head`:

<Execute command="head bla" />

To view the first 3 lines, use the `-n` parameter:

<Execute command="head -n 3 bla" />

To view the **last** 3 lines, we can instead use the `tail` command:

* <Execute command="tail -n 3 bla" />

Finally, the `grep` utility is a great way to filter lines in a file by patterns of interest:

* <Execute command="grep 'Promoter' bla" />

* <Execute command="grep -v 'Promoter' bla" />

* <Execute command="grep 'promoter' bla" />

* <Execute command="grep -i 'promoter' bla" />

* <Execute command="grep -i 'promoter' bla | wc -l" />

---

* Piping

---

* Output to a file with `>`

---

* <Execute command="abc=123" />

* <Execute command="echo $abc" />

* <Execute command="env" />

* <Execute command="USER=yourNameGoesHere" />

* <Execute command='echo "Hello $USER!"' />

* <Execute command="unset abc" />
