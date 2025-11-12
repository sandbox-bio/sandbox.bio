<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

Note that we may actually be interested in searching for the pattern in *only one* sequence of a record.
By default, `bqtools grep` will search for the patterns provided in *both* sequences of a record and return the record if it matches the pattern in *either* sequence.

We can limit our search to just the primary (R1) sequence by providing a pattern to the `-r/--reg1` flag or limit to just the extended (R2) sequence by providing a pattern to the `-R/--reg2` flag.

Let's test this out and find records that have the pattern `GGGG` in the primary sequence and `TTTTACGT` in the extended sequence.

<Execute command="bqtools grep merged.vbq -r GGGG -R TTTTACGT" />

Well that's not very many! How many exactly? Well lets count them!

<Execute command="bqtools grep merged.vbq -r GGGG -R TTTTACGT -C" />
