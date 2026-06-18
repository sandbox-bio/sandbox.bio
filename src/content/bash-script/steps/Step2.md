<script>
import Note from "$components/Note.svelte";
import Execute from "$components/Execute.svelte";
</script>

Although the script ran successfully, it's missing something _really_ important: the so-called "shebang" line, which looks like this:

```bash
#!/bin/bash
```

This line should be the first line of any Bash script, so the computer knows that it is in fact a Bash script. We got away with not specifying it because the default is Bash on sandbox.bio, but on many systems, the default is often `sh`, a shell that predates `bash`. While many things that work in `bash` work in `sh`, many other things will break on you in unexpected ways.

Next, **add that line at the very top of your script** and save your changes (Cmd + S on Mac, Ctrl + S on Windows).

<Note type="more">

You can put anything after the `#!`. For example, if you write `#!/usr/bin/python3` at the top of your Bash script, Python will try to run your script... and crash spectacularly of course because Python and Bash are different languages!

Another fun one is `#!/usr/local/bin/vim`, which will cause `./hello.sh` to open the file inside `vim`!

If you tried this and are now stuck in `vim`, I'm so sorry. Press `Escape`, `:`, `q`, and `Enter` until you find your way out.

</Note>
