<script> 
import Quiz from "$components/Quiz.svelte"; 
</script>

The `tree` command displays the tree-like organization of files and sub-directories contained in a particular directory.

In the example below, the `tree` command displays the content of the `/shared` directory limited only to directories (option `-d`) and with only two levels of sub-directories (option `-L 2`):

```bash
tree -d -L 2 /shared
```

These directories were created to store genome files of different species.

From the previous command we deduce the path from the root `/` to the `homo_sapiens` directory. This path is:

```bash
/shared/data/homo_sapiens
```

As stated previously, this path that starts with an `/` is an absolute path. Starting from the root `/`, we go through the `shared` then `bank` directories to reach the target `homo_sapiens` directory.

Paths are used in many Unix commands, such as the `ls` (that stands for **l**i**s**t) command:

```bash
ls /shared/data/homo_sapiens
```

This `ls` command lists the content of the specified directory (also named **argument** of the ls command).

<Quiz id="q1" choices={[
{ valid: false, value: "hg18"},
{ valid: true, value: "hg19"},
{ valid: false, value: "hg37"},
{ valid: true, value: "hg38"},
]}>
<span slot="prompt">
What does the command `ls /shared/data/homo_sapiens` return?
</span>
</Quiz>

Remark: Usually `hg` stands for **h**uman **g**enome and the number denotes the sequence version.
