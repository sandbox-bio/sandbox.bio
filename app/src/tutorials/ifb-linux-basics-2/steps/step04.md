<script>
import Quiz from "components/Quiz.svelte";
</script>

The `cd` command (**c**hange **d**irectory) can be used to visit another directory in the file tree. 
The target directory must be specified using an absolute or a relative path. 

To experiment with the `cd` command, run the following commands:

```bash
ls
cd /shared/data/nr
ls
cd ../homo_sapiens
ls 
```

<Quiz id="step04_q1" choices={[
	{ valid: false, value: "/shared/data/"},
	{ valid: false, value: "/shared/data/homo_sapiens"},
	{ valid: false, value: "../data/nr"},
	{ valid: false, value: "/shared/data/nr/homo_sapiens"},
	{ valid: true, value: "/shared/data/homo_sapiens"},
]}>
	<span slot="prompt">
		Could you guess the absolute path of your current working directory?
	</span>
</Quiz>

Verify the answer with `pwd`.

## Automatic completion

To go from your current working directory to a target directory, you must specify names of all intermediate directories. This can be time-consuming if the target directory is far away from your current directory. 
The key <kbd>Tab</kbd> triggers auto-completion. It means you just need to type the first letters of a directory, then <kbd>Tab</kbd>, to get its full name. If there is more than one file or directory starting with the same letter, auto-completion will complete the name as far as it can. If you type a second time <kbd>Tab</kbd>, auto-completion will show you the available options.

The <kbd>Tab</kbd> key is perhaps the most used key in Unix!

Use the <kbd>TAB</kbd> key and `cd` to go into the `/shared/data/bos_taurus/UMD3.1/star-2.7.2b/` directory.

<Quiz id="step04_q2" choices={[
	{ valid: false, value: "2"},
	{ valid: false, value: "3"},
	{ valid: false, value: "4"},
	{ valid: true, value: "5"},
]}>
	<span slot="prompt">
		How many files are in the `/shared/data/bos_taurus/UMD3.1/star-2.7.2b/` directory?
	</span>
</Quiz>
