<script>
import Quiz from "$components/Quiz.svelte";
</script>

## How to select part of a file

### _grep_

The `grep` command extracts lines matching a given pattern.
A pattern can be a simple word or a more general expression, often termed **regular expression** (see [here](https://librarycarpentry.org/lc-data-intro/01-regular-expressions/) to learn more about these).
For instance:

```bash
grep gene-SAOUHSC_00079 SAOUHSC.bed
```

<Quiz id="q1" choices={[
{ valid: false, value: "74750"},
{ valid: false, value: "94950"},
{ valid: true, value: "1561"},
{ valid: false, value: "750"},
]}>
<span slot="prompt">
Print the line that contains the gene-CD630_RS00010 gene name in the NC_009089.bed file. What is the starting position (given in the 2nd column) ?
</span>
</Quiz>

To count the number of lines containing the _gene_ word, just add the `-c` option to the `grep` command:

```bash
grep -c gene SAOUHSC.bed
```

<Quiz id="q2" choices={[
{ valid: false, value: "70"},
{ valid: false, value: "71"},
{ valid: true, value: "72"},
{ valid: false, value: "73"},
]}>
<span slot="prompt">
Count the number of lines containing the <i>cds</i> word in SAOUHSC.bed
</span>
</Quiz>

Here are other useful `grep` options :

- `-i`: searches the pattern in a case **i**nsensitive way
- `-n`: adds the line **n**umber at the beginning of the output line
- `-v`: prints the lines not containing the pattern (re**v**erse selection)

There are many others: try `man grep`
