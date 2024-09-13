<script>
import { Icon } from "sveltestrap";
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

Say you're studying the microRNA `miR-34a` and you want to extract its sequence in humans, you can use SeqKit's `grep` command:

<Execute command="seqkit grep --pattern 'hsa-mir-34a' hairpins.fa" />

To load sequences in bulk, SeqKit can read sequence IDs from a text file. For example, the file `ids.txt` contains the IDs of miR-34a across 5 species:

<Execute command="cat ids.txt" />

To extract just those sequences from `hairpins.fa`, use the `-f` flag:

<Execute command="seqkit grep -f ids.txt hairpins.fa" />

<Alert color="primary">
    <Icon name="lightbulb-fill" /> This is a useful tool because the alternative ways to achieve this on the command line would be more challenging. You could try using good old `grep` but that only matches the `>sequence_name` line. You could use `grep -A2` to also retrieve the sequence, but that assumes all sequences have 2 lines, which they do not.
</Alert>
