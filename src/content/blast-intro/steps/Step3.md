<script>
import Link from "$components/Link.svelte";
import Image from "$components/Image.svelte";
import Execute from "$components/Execute.svelte";
</script>

To put these various tools and options to use, let's consider using `blastp` to look for proteins that are similar in sequence to other proteins in the yeast exome. For this tutorial, we'll use the FASTA file `orf_trans.fasta`:

<Execute command="head orf_trans.fasta" />

In order to find sequences that are similar to others, we're going to want to `blastp` this file against itself. So, we'll start by creating a database of these sequences:

<Execute command={`makeblastdb \\ -in orf_trans.fasta \\ -out orf_trans \\ -dbtype prot \\ -title "Yeast Open Reading Frames" \\ -parse_seqids`} />

This will take a few seconds to complete.

**Note:** The `-parse_seqids` flag indicates that the sequence IDs from the FASTA file should be included in the database so that they can be used in outputs as well as by other tools

Once you run the command above, notice the new database files that were created:

<Execute command="ls -l" />
