<script>
import Quiz from "$components/Quiz.svelte";
import Image from "$components/Image.svelte";
</script>

Using the filesystem diagram below,
if `pwd` displays `/Users/backup`,
and `-r` tells `ls` to display things in reverse order,
what command(s) will result in the following output:

```bash
pnas_sub/ pnas_final/ original/
```

<Image src="https://swcarpentry.github.io/shell-novice/fig/filesystem-challenge.svg" alt="File system diagram" />

<Quiz id="q2.3" choices={[
{ valid: false, value: "ls pwd"},
{ valid: true, value: "ls -r -F"},
{ valid: true, value: "ls -r -F /Users/backup"},
]}>
<span slot="prompt">
</span>
</Quiz>
