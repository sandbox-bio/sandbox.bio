<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

To use <Link href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-versus-local-alignment">local alignment</Link> to align some longer reads included with Bowtie 2, stay in the same directory and run:

<Execute command="bowtie2 \ --local \ -x $REF \ -U longreads.fq \ -S eg3.sam" />

This aligns the long reads to the reference genome using local alignment, with results written to the file `eg3.sam`.
