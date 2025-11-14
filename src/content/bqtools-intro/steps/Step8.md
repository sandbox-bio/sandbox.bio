<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

One of the most useful things about BINSEQ files is that cheap per-sequence tasks can be trivially parallelized.
This is great, because many problems in bioinformatics belong to a class of problems that are <Link href="https://en.wikipedia.org/wiki/Embarrassingly_parallel">embarrassingly parallel</Link>.

One such task is searching for specific subsequences within each sequence.
This is equivalent to performing a `grep` command on each sequence.

This is such a common operation that `bqtools` provides a built-in command for it: `bqtools grep`.

Let's look for the subsequence `ACGTACGT` in our records.

<Execute command="bqtools grep merged.vbq ACGTACGT" />

Looks like we got some hits! You'll see that `bqtools` will highlight the matches for you in the sequence.

Let's actually just count the number of records we got with this search using the `-C/--count` flag.

<Execute command="bqtools grep merged.vbq ACGTACGT -C" />

We can also add additional patterns to search for in our records.
By default the logic for accepting a read with multiple patterns follows an `AND` logic, which means that a read must contain all patterns to be accepted.

Let's find all the records that contain _both_ `ACGTACGT` and `TATA` and output as a fastq using the `-f/--format` option:

<Execute command="bqtools grep merged.vbq ACGTACGT TATA -fq" />

> Note: to search for patterns using `OR` logic, you can use the `--or-logic` flag.

`bqtools grep` is a fully <Link href="https://en.wikipedia.org/wiki/Regular_expression">Regex</Link> compliant tool, which means you can use regular expressions to search for sequences.

Let's take a look at a more complicated query - perhaps a sequence that begins with `ACGT` but then has either `TA` or `GC`, and finally ends with `CCGG`.

<Execute command='bqtools grep merged.vbq "ACGT(TA|GC)CCGG"' />
