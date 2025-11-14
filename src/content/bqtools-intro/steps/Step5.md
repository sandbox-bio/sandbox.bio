<script>
import Execute from "$components/Execute.svelte";
import Link from "$components/Link.svelte";
</script>

`bqtools` has some nice features that make it easy to work with existing sequencing data.
One of these is recursive encoding - which will convert all the sample pairs in a directory (and subdirectories) to their BINSEQ representation.

Let's try it out with all the sample pairs we have in the `data` directory:

<Execute command="bqtools encode --recursive fastq" />

Oh, there is a warning!

Let's take a look at the files that were created:

<Execute command="ls -lh fastq/*.vbq" />

We can see that we have 8 VBQ files with an R1 or R2 - but we've learned we do not need to keep these files separate.

The `--recursive` flag tells `bqtools` to convert every *individual* file to a BINSEQ file.
But it also correctly identified that the files in this directory appear to be paired.
So it converted each individual file and did not pair them - but it warned us that this might not be what we wanted.

Let's remove those individual files that were not correctly paired and try again:

<Execute command="rm -v fastq/*.vbq" />

Now let's try again with the `--paired` flag:

<Execute command="bqtools encode --recursive --paired fastq" />

Great! Now we have exactly what we want - 4 files that are paired:

<Execute command="ls -lh fastq/*.vbq" />

> Note: the `--recursive` flag will convert all files in a directory and its subdirectories to BINSEQ files. If you want to limit the depth of the recursion you can use the `--depth` flag.
