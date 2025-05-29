<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

In addition to the reference-free alignments, SKA can also be used to produce alignments against a reference genome. As an example, let's use a reference genome for _E. coli_, it is already downloaded for you, but you can download it in your computer from <Link href="https://www.ebi.ac.uk/ena/browser/api/fasta/U00096.3?download=true">this ENA entry</Link> if you want.

To do so, we need to create a new alignment, but adding the reference to it. We have provided another input file list, `assemblies/ska_input_ref.tsv`, that replaces one of the previous entries with the reference, so let's run again the alignment with:

<Execute command={"ska build \ -f assemblies/ska_input_ref.tsv \ -k 31 \ -o output/reference_index"} />

And, with it, we can align against the reference genome using `ska map`:

<Execute command={"ska map \ -o output/reference_map.aln \ --ambig-mask assemblies/GCA_000005845.2.fna.gz \ output/reference_index.skf"} />

Compared to `ska align`, `ska map` keeps all bases in the reference sequence and replaces the sites it cannot find in the indexed assemblies with gaps. We can create the tree again with our script:

<Execute command="python3 create_tree.py output/reference_map.aln" />

Here we can see that the reference is not closely related to any of the other entries. This makes complete sense, as the reference we took comes from the _E. coli_ K-12 strain, isolated one century ago on the other side of the world (Palo Alto, California, USA)!
