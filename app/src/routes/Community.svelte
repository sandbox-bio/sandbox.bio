<script>
import Listings from "./components/Listings.svelte";
import { tutorials } from "./stores/tutorials";
</script>

<h3 class="my-3">Community Tutorials</h3>

<h5 class="bg-primary bg-opacity-75 p-2 text-white">French Institute of Bioinformatics</h5>
<img src="/logo-ifb.svg" alt="Logo of the IFB" width="300" />

<Listings title={null} items={$tutorials.filter(t => t.id.startsWith("ifb-"))} showUnlisted={true} />
