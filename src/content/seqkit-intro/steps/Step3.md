<script>
import { Icon } from "sveltestrap";
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
</script>

The `seqkit seq` command is used to extract, filter and format your FASTA and FASTQ files.

For example, to extract the sequence names from a FASTA file:

<Execute command="seqkit seq --name hairpins.fa | head" />

If your FASTA is formatted such that the sequence name contains an ID followed by a space and more information, then you can extract just those IDs using `--only-id`:

<Execute command="seqkit seq --name --only-id hairpins.fa | head" />

<hr />

If you are interested in only analyzing hairpins that are >300bp long, use the `--min-len` to filter out shorter sequences:

<Execute command="seqkit seq --min-len 300 hairpins.fa  | seqkit stats" />

You can also filter out long sequences with `--max-len`, and for FASTQ files, filter out reads with a certain average quality with `--min-qual` and `--max-qual`.

<Alert color="primary">
    <Icon name="question-circle-fill" /> How would you convert the RNA sequences in `hairpins.fa` to DNA using SeqKit? Use the <Link href="https://bioinf.shenwei.me/seqkit/usage/#seq">manual</Link> as a reference.
</Alert>
