<script>
import Listings from "$components/Listings.svelte";
import { tutorials } from "$stores/tutorials";
</script>

<h3 class="my-3">Community Tutorials</h3>

<h5 class="bg-primary bg-opacity-75 p-2">
	<a href="https://www.france-bioinformatique.fr/en/home/" class="text-white text-decoration-none" target="_blank" rel="noreferrer">
		French Institute of Bioinformatics
	</a>
</h5>

<a href="https://www.france-bioinformatique.fr/en/home/" target="_blank" rel="noreferrer">
	<img src="/logo-ifb.svg" alt="Logo of the IFB" width="250" />
</a>

<Listings title={null} items={$tutorials.filter((t) => t.id.startsWith("ifb-"))} showUnlisted={true} />

<p class="text-muted small mt-3">
	IFB tutorials are licensed under a
	<a href="https://creativecommons.org/licenses/by-sa/4.0/legalcode" target="_blank" rel="noreferrer"> CC-BY-SA </a>
	license.
</p>
