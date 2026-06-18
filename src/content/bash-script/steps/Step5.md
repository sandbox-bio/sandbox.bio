<script>
import Note from "$components/Note.svelte";
import Execute from "$components/Execute.svelte";
</script>

Currently, the script doesn't give an error if you only give it one parameter:

<Execute command="./hello.sh Hello" />

Let's change that by using the handy `?` notation when defining variables:

```bash|editor
#!/bin/bash
greeting=${1?Error: please specify a greeting}
name=${2?Error: please specify a name}

echo "$greeting $name"
```

Now try again:

<Execute command="./hello.sh Hello" />

This works, but a pattern you'll see in most command-line programs is to output usage information when there's any error in the way the script was called:

```bash|editor
#!/bin/bash
usage="Usage: ./hello.sh greeting name"
greeting=${1?$usage}
name=${2?$usage}

echo "$greeting $name"
```

<Execute command="./hello.sh Hello" />

This is useful because you only need to write the usage once, and not one error per parameter. It's also good because it gives others information about other required parameters, so users of the script don't need to fix their inputs, run the code, fix their inputs again, etc.
