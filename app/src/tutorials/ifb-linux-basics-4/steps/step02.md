<script>
import Quiz from "components/Quiz.svelte";
import Execute from "components/Execute.svelte";
</script>

The standard command input is named **stdin**.

<img src="/data/linux_basics_session04/stream_in_out.png" style="max-width:100%" alt="input stream of a command">

By default, **stdin** is set to the keyboard. But you can change this behavior and read **stdin** from a file. To do this, you need to use the `<` operator.

The `tr` command translates (or deletes) characters from a text supplied as input. If we want to convert the `NC_009089.fasta` file that contains a DNA sequence into an RNA sequence (replacing T with U) and switch from upper to lower cases, we can use the following command:

<Execute command="tr AGCT agcu < NC_009089.fasta " />

Remark: `cat toto` is equivalent to `cat < toto` but `cat < toto` is seldom used.

<img src="/data/linux_basics_session04/stream_infile_out.png" style="max-width:100%" alt="input stream of a command">
