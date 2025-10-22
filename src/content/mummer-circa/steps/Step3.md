<script>
    import Execute from "$components/Execute.svelte";
</script>

In Circa, we need to give it the names and sizes of the chromosomes to show.

When showing a single genome, the .fai index file is enough, but when showing two genomes, we need to create a reference file that includes the sizes and orientations of the chromosomes in the order we want to show them.

First, let's index our FASTA files:

<Execute command={`samtools faidx $REF`} />
<Execute command={`samtools faidx $QUERY`} />

The simplest reference file is just the .fai files concatenated together:

<Execute command={`cat $REF.fai $QUERY.fai > $NAME.ref.tsv`} />
The result looks like this, with the first two columns being the names and sizes.
<Execute command="cat $NAME.ref.tsv" />

That's all we need to set up the coordinates for visualization.

Next, we will use this reference file and the alignments from the previous step to make a circos plot in Circa.
