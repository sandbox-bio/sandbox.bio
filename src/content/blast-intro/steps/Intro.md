<script>
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import Image from "$components/Image.svelte";
</script>

<Alert>
This tutorial is an interactive version of the <Link href="https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/">BLAST chapter</Link> of the open-access textbook **A Primer for Computational Biology**. The book is licensed under a <Link href="https://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons license</Link>.
</Alert>

BLAST, or Basic Local Alignment Search Tool, is a nearly ubiquitous tool for sequence alignment.

Given one or more _query_ sequences (usually in FASTA format), BLAST looks for matching sequence regions between them and a _subject_ set:

<Image src="/data/blast-intro/blast.png" alt="How blast works: look for matching sequence between sets of queries and subjects." />

A sufficiently close match between subsequences (denoted by arrows in the figure above, though matches are usually longer than illustrated here) is called a **high-scoring pair (HSP)**, while a query sequence is said to hit a target sequence if they share one or more HSPs. Sometimes, however, the term "hit" is used loosely, without differentiating between the two.

Each HSP is associated with a **bitscore** that is based on the similarity of the subsequences as determined by a particular set of rules. Because in larger subject sets some good matches are likely to be found by chance, each HSP is also associated with an **E value**, representing the expected number of matches one might find by chance in a subject set of that size with that score or better. For example, an E value of 0.05 means that we can expect a match by chance in 1 in 20 similar searches, whereas an E value of 2.0 means we can expect 2 matches by chance for each similar search.
