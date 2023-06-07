<script>
import Quiz from "components/Quiz.svelte";
</script>

## Move (or rename) files and directories

The `mv` (**m**o**v**e) command takes 2 paths as arguments:

```bash
mv <source_path> <destination_path>
```

It moves the **source** to the **destination**.  
It works for files or directories.  
It is also used to rename files or directories.

Try this `mv` command to rename the file from `second_file.txt` to `2nd_file.txt` and to move it towards your HOME directory: 

```bash
cd
tree
mv test/second_file.txt 2nd_file.txt
tree
```

## Delete files and directories

The `rm` (**r**e**m**ove) command deletes files or directories.

Use `rm` to delete the file named `second_file.txt` from the directory `~/test`. Also use `tree` to check the organisation of files and directory from the working directory:

```bash
cd
tree
rm test/first_file.txt
tree
```

To delete a directory, you need to use the `rm` command with the option `-r`:

```bash
rm -r <path_to_a_directory_to_delete>
```

Be very careful with this `rm` command. There is no way to recover your deleted files in Unix!

<Quiz id="step07_q1" choices={[
    { valid: false, value: "`mv` applies to files or directories while `rm` applies to directories only"},
	{ valid: true, value: "`mv` requires 2 paths while `rm` requires only one path"},
]}>
	<span slot="prompt">
		The `mv` command differs from the `rm` command by (select the right proposal):
	</span>
</Quiz>
