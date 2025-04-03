<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

The terminal is a text-based interface that interprets commands and outputs the result to your screen.

The first command we'll try is `echo`, which simply returns the string you provide it.

Click the box below to execute the command (or type it in manually to improve the retention of what you learn):

<Execute command='echo "Hello World"' />

Try variations of that command with your own strings.

Just like the graphical interface you use for exploring files on your computer, you always work within one folder when you're using the terminal.

For most of the tutorials on sandbox.bio, you won't have to worry about which folder you're in, but you can use `pwd` to **p**rint the **w**orking **d**irectory:

<Execute command="pwd" />

To list the files inside this folder, use `ls`:

<Execute command="ls" />

You should see three files, one of which is `orders.tsv`, which contains data about restaurant <Link href="https://github.com/TheUpshot/chipotle/">take-out orders</Link>.

In the next step, we'll start exploring that data.
