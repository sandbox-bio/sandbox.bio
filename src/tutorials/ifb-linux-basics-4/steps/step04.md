<script>
import Execute from "components/Execute.svelte";
</script>

As we've seen so far, a file can be redirected to the **stdin** of a command using the `<` operator, while the stdout of a command can be redirected to a file using the `>` operator (or `>>` to add lines).

We want to count the non-redundant list of entries in the fourth column of the SAOUHSC.bed file. 

Use the **cut** command to extract the 4th column from the bed file SAOUHSC.bed and create a file named `SAOUHSC_c4.bed`. The column to be cut must be specified using `-f 4` or `--fields 4` if using the long form of the argument. Note that the file is supposed to be tabulated by default (columns must be separated by a '\t'). This behavior can be modified using the `-d/--delimiter` argument.  

<Execute command="cut -f 4 SAOUHSC.bed > SAOUHSC_c4.bed" />

We can search for the list of entries in column 4 using the sort command (which by default performs an alphanumeric sort) combined with the `-u` (or --unique in its long form) to ensure that the list is not redundant.


<Execute command="sort -u SAOUHSC_c4.bed > SAOUHSC_c4_uniq.bed" />

If we want to count the number of non-redundant entries in the SAOUHSC_c4_uniq.bed file, we now need to use the `wc` command. By default, the `wc` command counts the number of lines, words and bytes in a file. We can use the `-l` or `--lines` arguments to simply count the number of lines.

<Execute command="wc -l SAOUHSC_c4_uniq.bed" />

As you may have noticed, this simple step required us to create several intermediate files that won't be particularly useful later on. In the next step, we'll see that we can get rid of these files by using the pipe `|` operator.
