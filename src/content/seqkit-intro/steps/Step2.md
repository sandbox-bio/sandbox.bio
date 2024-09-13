<script>
import { Icon } from "sveltestrap";
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

If you're testing a new tool or writing a new algorithm, working with large data files can slow you down because long runtimes make it difficult to quickly iterate on your work.

This is where downsampling (or subsampling) comes in. In the previous step, we saw that `hairpins.fa` has ~3.1K sequences. Let's sample 10% of that file and save it as its own FASTA file:

<Execute command="seqkit sample --proportion 0.1 hairpins.fa > sampled.fa" />

Instead of a fraction, to obtain a number of sequences sampled, use the `--number` flag:

<Execute command="seqkit sample --number 10 hairpins.fa > sampled.fa" />

<Alert color="warning">
    <Icon name="lightbulb-fill" /> Depending on the random seed, you may not always obtain exactly the number of sequences requested. For example:<br /><br />

<Execute command="seqkit sample --number 10 --rand-seed 123 hairpins.fa > sampled.fa" />

See the <Link href="https://bioinf.shenwei.me/seqkit/note/#effect-of-random-seed-on-results-of-seqkit-sample">SeqKit manual</Link> for details.

</Alert>
