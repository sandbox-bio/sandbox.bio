<script>
import Note from "$components/Note.svelte";
import Execute from "$components/Execute.svelte";
</script>

Scripts are most useful when you can change their behavior without modifying the script itself, which is what command-line parameters are for.

For example, let's update our script to take 2 parameters: a greeting and a name, and have it output the greeting followed by the name.

We can do that by storing the value of each command-line parameter in its own variable. The first parameter is represented by the variable `$1`, and the second parameter by `$2`:

```bash|editor
#!/bin/bash
greeting=$1
name=$2

echo "$greeting $name"
```

Try running this script:

<Execute command="./hello.sh 'Good morning' Robert" />

<Note type="whatif">

What happens if you don't put quotes around `Good morning` when calling the script?

</Note>

<Note type="more">

In Bash, you can define variables using <Execute inline command="myname=John" />, then use those variables using `$myvalue`, e.g. <Execute inline command="echo $myname" />

Make sure there are no spaces around the equal sign! For example, `myname = value` in Bash would mean: run the program `myname` and give it two parameters: `=` as the first parameter, and `value` as the second!

Inside a script, the variables `$1`, `$2`, ... represent the `n`th command-line argument, and `$0` is the name of the program you ran (try adding `echo $0` to your script).

</Note>
