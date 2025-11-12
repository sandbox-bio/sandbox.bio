<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

A feature that doesn't exist on most grep-like tools is the ability to count the number of matches found *per-pattern*.
This is actually quite a useful feature for bioinformatics and shows up often when exploring data.

Let's pretend that we have a collection of patterns - perhaps these are sections of a transcript, or a list of cell-barcodes, or perhaps sgRNA protospacers - and we want to count how many times each of them appear in our dataset.

Let's take a look at these patterns:

<Execute command="cat patterns/patterns.txt" />

`bqtools grep` will accept a text file containing a list of patterns to search for - so let's count how many times each of these patterns appear in our dataset:

<Execute command="bqtools grep merged.vbq --file patterns/patterns.txt -P" />

> Note that we use the `-P/--pattern-count` flag to enter the pattern-counting mode of `bqtools grep`.

This will give us a table with 3 columns that reflect the pattern, the number of sequence matches, and the fraction of all records that match that pattern.

> If you have a large number of patterns and they are not regular expressions, try using the `-x/--fixed` flag which makes use of an optimized algorithm to count the number of matches.
