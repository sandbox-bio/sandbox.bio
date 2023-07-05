<script>
import Link from "$components/Link.svelte";
import Alert from "$components/Alert.svelte";
import Execute from "$components/Execute.svelte";
import Listings from "$components/Listings.svelte";
import { tutorials } from "$stores/tutorials";
</script>

This concludes our introduction to the command line!

You are now equipped to tackle the other tutorials on sandbox.bio.

<a href="/tutorials?id=awk-intro" class="btn btn-primary px-4 me-md-2 fw-bold">Keep going &rarr;</a>

<Listings title="More Tutorials" colXxl={12} colMd={12} colLg={12} items={$tutorials} skip={["terminal-basics"]} />
