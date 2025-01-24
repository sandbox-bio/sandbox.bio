<script>
import Link from "$components/Link.svelte";
import Execute from "$components/Execute.svelte";
</script>

**1. The sequencing reads**

We will be using 10 HIV-1 whole genome sequences collected from samples from real people! These sequences can be found in the file `hiv1_sequences.fas`.

Try <Execute command="less hiv1_sequences.fas" inline /> to
take a look at the HIV-1 sequences we will be using.

**2. Multiple sequence alignment (MSA)**

We will begin by performing multiple sequence alignment for our HIV sequences. Use <Execute command="mafft hiv1_sequences.fas > hiv1_sequences.MSA.fas" inline /> to run MAFFT to perform MSA.
