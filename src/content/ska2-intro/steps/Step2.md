<script>
import Execute from "$components/Execute.svelte";
</script>

The SNP alignment can be extracted from the .skf file produced in the previous step by running

<Execute command={"ska align \ --min-freq 1 \ --filter no-filter \ -o output/ska_alignment.aln \ output/ska_index.skf"} />

which will again take a bit (~2 min) to run. In the above command, the `--min-freq` option sets the maximum number of missing sites in the alignment to 1 and the `--filter no-filter` option specifies that no extra filtering should be performed.

Other options include filtering all constant sites with `--filter no-const` or all constant and ambiguous sites with `--filter no-ambig-or-const`. Using the last option is equivalent to treating any ambiguous site as an N.
