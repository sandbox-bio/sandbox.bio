<script>
import Execute from "$components/Execute.svelte";
import Note from "$components/Note.svelte";
</script>

Let's start out with a "simple" example. In the code editor on your right, you'll find the script `hello.sh`, which outputs `Hello world` to the screen.

The more detail-oriented among you may have noticed the quotes around the word "simple", and you'd be right to think something's up because this script doesn't run!

Try it yourself:

<Execute command="./hello.sh" />

You'll see a `Permission denied` error, but don't worry, this is an error most of us run into, even after years of using the command line.

It happens because when you create a new file in Linux, you only have permissions to read (`r`) and write (`w`) it, but not execute (`x`) it.

For that, you need to explicitly grant yourself (`u`) permission to do so:

<Execute command="chmod u+x hello.sh" />

Then run the script again:

<Execute command="./hello.sh" />

Now it's working as expected!

<Note>

Most people just use `chmod +x` without the `u`, but keep it mind that gives `execute` permissions to all users on the system, which is fine on your laptop but can be a problem on shared servers, where you might not want anyone else to run your scripts.
</Note>

<Note type="more">

`chmod` stands for <span class="text-primary fw-bold">ch</span>ange <span class="text-primary fw-bold">mod</span>e, and <Execute inline command="ls -l" /> shows the mode of each file and folder in the very first column, as `r`, `w` and `x` flags.

That string in the first column (`-rwxr--r--`) actually breaks down into 4 parts as `-`, `rwx`, `r--`, and `r--` as follows:

- A character that is set to `d` if the item is a directory
- 3 characters that show which permissions you have (in this case read, write, execute)
- 3 characters that show permissions of users in your group (in this case, only read access)
- 3 characters that show permissions of all other users (in this case, only read access)

This will look more interesting if you try it on your work servers or university cluster where there are more users and groups set up.

</Note>
