<script>
import Quiz from "components/Quiz.svelte";
import Execute from "components/Execute.svelte";
</script>

## the standard output strem

By default, as you've seen so far, the result of a Unix command is printed on the screen.

As an example, we will extract genomic locations related to *gene-SAOUHSC_00079* entry from the `SAOUHSC.bed` file with the `grep` command:
First, check that you have access to the `SAOUHSC.bed` file using the `ls` command, and then extract the location with the following instructions:

<Execute command="ls" />

<Execute command="ls Data" />

<Execute command="cd Data" />

<Execute command="grep SAOUHSC_00079 SAOUHSC.bed" />

The result of the `grep` command is displayed on the terminal.

Here are some vocabulary definitions:

📕 the result of a command is also called **output** 

📕 the **st**andar**d** **out**put of a command is named **stdout**.

The following diagram illustrates the output stream of a command:

<img src="/data/linux_basics_session04/stream_out.png" style="max-width:100%" alt="stream_out">

By default, **stdout** is set to the screen.

## changing the standard output stream

You can modify this behavior and print **stdout** to a file.
To do so, you need to use the `1>` that can be abbreviated to `>`:

<Execute command="grep gene-SAOUHSC_00079 SAOUHSC.bed > gene.bed" />

Look, you've created a new file named `gene.bed`

<Execute command="ls" />

You can view its contents using the `cat` command:

<Execute command="cat gene.bed" />

The content of this new file is identical to the result of the `grep` command.

The `>` symbol is one of the **redirection** operators.

The following figure illustrates the **stdout** redirection to a file:

<img src="/data/linux_basics_session04/stream_outfile.png" style="max-width:100%" alt="stream_outfile">

⚠️ if the file already exists, it’s contents will be replaced by the output of your instruction.

If you execute the same `grep` instruction as before but search for a different gene, the output file will be overwritten:

<Execute command="grep gene-SAOUHSC_00078 SAOUHSC.bed > gene.bed" />

<Execute command="cat gene.bed" />

If you want to store both gene locations in a single file, you may use the `>>` operator, which appends the output of your command to the end of an existing file.

<Execute command="grep gene-SAOUHSC_00079 SAOUHSC.bed > gene.bed" />

<Execute command="grep gene-SAOUHSC_00078 SAOUHSC.bed >> gene.bed" />

<Execute command="cat gene.bed" />

<Quiz id="q1" choices={[
         { valid: false, value: "grep foo file01 > file02"},
         { valid: true, value: "grep foo file01 >> file02"},
         { valid: false, value: "grep foo file02 >> file01"},
	 { valid: false, value: "grep file01 file02 > foo"},
]}>
        <span slot="prompt">
	Which command appends its result at the end of the file02 file ?
        </span>
</Quiz>
