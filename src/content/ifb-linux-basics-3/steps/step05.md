<script>
import Quiz from "$components/Quiz.svelte";
import Execute from "$components/Execute.svelte";
</script>

## a bed file

The `bed` file format is widely used in bioinformatics.
It can be used to represent known sequence information, by their position on chromosomes.
For example, see the content of the `NC_009089.bed` file with the `cat` command.

<Execute command="cat NC_009089.bed" />

In the 6-column bed file `NC_009089.bed`, we find in order: the sequence identifier (here chromosome NC_009089, version 1 in the column 1), the start (column 2) and end (column 3) positions of the annotation, the annotation identifier (column 4, with genes and coding sequences, cds), column 5 is sometimes used to indicate a confidence score (here there is a dot meaning that this column is not used) and, in column 6, the strand (forward `+` or reverse `-` direction).

## _wc_

`wc` is the **w**ord **c**ount command.
It counts the number of lines, words, and characters in files:

<Execute command="wc NC_009089.bed" />

The _NC_009089.bed_ has from left to right: 8 lines, 48 words and 364 characters.

To output only the number of lines, run `wc` with the _-l_ option :

<Execute command="wc -l NC_009089.bed" />

<Quiz id="q1" choices={[
{ valid: false, value: "5"},
{ valid: false, value: "28"},
{ valid: false, value: "67"},
{ valid: true, value: "150"},
]}>
<span slot="prompt">
How many lines has the SAOUHSC.bed file?
</span>
</Quiz>
