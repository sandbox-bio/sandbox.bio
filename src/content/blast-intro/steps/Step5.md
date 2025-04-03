<script>
import Link from "$components/Link.svelte";
import Image from "$components/Image.svelte";
import Execute from "$components/Execute.svelte";
</script>

We can see with `tail` that the output file contains the columns we specified interspersed with the comment lines provided by `-outfmt 7`:

<Execute command={`tail yeast_blastp_yeast_top2.txt`} />

In the output above, `YHR007C` has an HSP with itself (naturally), but also one with `YDR402C`.
