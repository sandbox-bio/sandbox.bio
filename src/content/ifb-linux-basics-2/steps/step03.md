<script>
import Quiz from "$components/Quiz.svelte";
</script>

Run the `pwd` command in the right panel.

The output is the absolute path `/root/tutorial` and means that you are currently located in `/root/tutorial`.

Remark: In a Unix system, the administrator (the boss) is called root. And you are presently in its personal directory that is also called `root`!

If you type the `ls` command, you should see a sud-directory called `test`:

```bash
ls
```

From the present current working directory, we would like to see what is inside a sub-directory called `test`.
To represent the current working directory, we need the symbol `.` (dot).
The relative path to the `test` sub-directory is `./test`.

```bash
ls ./test
```

By default, relative paths start from the current working directory, so `./` could be omitted:

```bash
ls test
```

The `..` operator is handy to write a path relative to a directory.
It means _one level up in the directory tree_.
For instance from the `/root/tutorial` directory where you are, you may list the content of the `/root` directory using:

```bash
ls ..
```

The same result would be obtained here using an absolute path:

```bash
ls /root/
```

In an other example, if you are located in the `/root` directory, you could list the content of `/tmp` with a relative path:

```bash
ls ../tmp
```

The same result would be obtained with the absolute path:

```bash
ls /tmp
```

<Quiz id="q1" choices={[ { valid: false, value: "/"},
{ valid: false, value: "/root"},
{ valid: true, value: "/root/tutorial"}, ]}>
<span slot="prompt">
If your current working directory is `/root/tutorial/homo_sapiens`, to which absolute path refers the path `..` ?
</span>
</Quiz>

<Quiz id="q2" choices={[ { valid: true, value: "absolute"},
{ valid: false, value: "relative"}, ]}>
<span slot="prompt">
Which type of path is `/root/tutorial`?
</span>
</Quiz>

<Quiz id="q3" choices={[ { valid: false, value: "homo_sapiens/hg19/fasta"},
{ valid: false, value: "../../hg19/fasta"},
{ valid: true, value: "../homo_sapiens/hg19/fasta"}, ]}>
<span slot="prompt">
If your current working directory is `/root/tutorial/bos_taurus` what is the relative path to `/root/tutorial/homo_sapiens/hg19/fasta`?
</span>
</Quiz>
