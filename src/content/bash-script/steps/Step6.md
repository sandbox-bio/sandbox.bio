<script>
import Note from "$components/Note.svelte";
import Execute from "$components/Execute.svelte";
</script>

Instead of making the `name` required, we can also make it optional and give it a default value using the weird `:-` notation (there's a lot of weird things about Bash, wait till you hear about arrays):

```bash|editor
#!/bin/bash
usage="Usage: ./hello.sh greeting [name]"
greeting=${1?$usage}
name=${2:-World}

echo "$greeting $name"
```

Now run the script again without the second parameter to see if it picks up the default value:

<Execute command="./hello.sh Hello" />

But if you specify a second parameter, it overwrites the default, as expected:

<Execute command="./hello.sh Hello Robert" />

<Note>

I won't cover it here, but it is possible for a Bash script to handle parameters that are not based on position, with flags like `./hello.sh -n Robert -g Hello`, where the order doesn't matter anymore. I won't show it here but you can use [getopts](https://stackoverflow.com/a/16496491) to do so.

I will say though, this gets really messy in Bash, so if you're finding yourself using named parameters in a Bash script, it might be time to switch to another language!

</Note>
