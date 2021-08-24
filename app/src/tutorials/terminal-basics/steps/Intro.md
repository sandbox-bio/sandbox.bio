<script>
import Link from "components/Link.svelte";
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

You should see the file `orders.tsv`, which contains data about restaurant <Link href="https://github.com/TheUpshot/chipotle/">take-out orders</Link> placed in 2015.

In the next step, we'll start exploring that data.

---

To view the first 10 lines of our orders dataset, use `head`:

<Execute command="head orders.tsv" />

To view the first 3 lines, use the `-n` parameter:

<Execute command="head -n 3 orders.tsv" />

To view the **last** 3 lines, we can instead use the `tail` command:

<Execute command="tail -n 3 orders.tsv" />

From that, we can infer that there is information about `1834` orders, and that the columns of this file contain the order ID, the quantity, and the item ordered.

If we're just interested in the lines in `orders.tsv` that contain burrito orders, we can use `grep` to filter lines using a pattern of interest:

<Execute command='grep "Burrito" orders.tsv' />

Note that `grep` is case-sensitive:

<Execute command='grep "burrito" orders.tsv' />

For convenience, we can ask `grep` to ignore case:

<Execute command='grep -i "burrito" orders.tsv' />

Another command pattern of inquiry is to find lines that do **not** match a pattern. To find all orders that aren't burritos:

<Execute command='grep -v "Burrito" orders.tsv' />

---

Next, let's cover **piping**, which helps you string together many command-line tools, such that the output of one is the input of the other.

For example, to find all orders that aren't burritos, and only display the last 3, we can pipe (i.e. `|`) the output of `grep` to `tail`:

<Execute command='grep -v "Burrito" orders.tsv | tail -n 3' />

We can also use the `wc -l` command to only **count** the number of lines without displaying them:

<Execute command='grep "Chicken Burrito" orders.tsv | wc -l' />

<Execute command='grep "Steak Burrito" orders.tsv | wc -l' />

---

We can store the result of our analyses in files using the `>` operator.

For example, to store all chicken burrito orders into `burritos.tsv`:

<Execute command='grep "Chicken Burrito" orders.tsv > burritos.tsv' />

Now notice that <Execute command="ls" inline /> shows the newly-created file!

<Alert color="warning">

Be careful when using `>`. If you have a FASTA file and you want to extract lines with the character `>`, make sure to put it in quotes! Otherwise, grep interprets `>` as an operator and will truncate your FASTA file.

This is a rite of passage in bioinformatics but can be avoided if you prepend `cat` to your pipelines, e.g. `cat ref.fa | grep ">"` (this is why there is no such thing as a _useless use of cat_).

<!-- Let's illustrate this pitfall in our sandbox. If we grep for `>` with quotes, it works as expected:
<Execute command='grep ">" copy.fa' />
Whereas not quoting results in a truncated file!
<Execute command='grep > copy.fa' />
You can verify that using `ls`:
<Execute command='ls' /> -->
</Alert>

---

Finally, let's explore environment variables. These are variables you can define in the terminal:

<Execute command="abc=123" />

Make sure there are no spaces surrounding the equal sign, otherwise the terminal treats the variable name as a command!

To display the content of a variable, use `echo`:

<Execute command="echo $abc" />

To delete a variable, use `unset`:

<Execute command="unset abc" />

To list all variables available in your environment:

<Execute command="env" />

Note that there's a variable called `USER` that we use while displaying the command-line prompt.

You can modify this variable to customize your environment:

<Execute command="USER=yourNameGoesHere" />

Now your prompt should update to reflect your name instead of `guest`!
