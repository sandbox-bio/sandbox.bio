<script>
import { Icon } from "sveltestrap";
import Link from "components/Link.svelte";
import Alert from "components/Alert.svelte";
import IGVUpdateBtn from "components/IGVUpdateBtn.svelte";
</script>

Navigate to a narrow window on chr21: `21:19,480,041-19,480,386`

<IGVUpdateBtn locus="21:19,480,041-19,480,386" />

You will see reads represented by bars stacked on top of each other, where they were aligned to the reference genome. The reads are pointed to indicate their orientation (i.e. the strand on which they are mapped).

Click on any read and notice that a lot of information is available. Once you select a read, you will learn what many of these metrics mean, and how to use them to assess the quality of your datasets.

At each base that the read sequence **mismatches** the reference, the colour of the base represents the letter that exists in the read (using the same colour legend used for displaying the reference).
