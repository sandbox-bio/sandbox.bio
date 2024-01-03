<script>
import TutorialList from "$src/components/TutorialList.svelte";
import { config as ifblinuxbasics1 } from "$content/ifb-linux-basics-1/config.js";
import { config as ifblinuxbasics2 } from "$content/ifb-linux-basics-2/config.js";
import { config as ifblinuxbasics3 } from "$content/ifb-linux-basics-3/config.js";
</script>

<svelte:head>
	<title>Community - sandbox.bio</title>
</svelte:head>

<h3 class="my-3">Community Tutorials</h3>

<h5 class="bg-primary bg-opacity-75 p-2">
	<a href="https://www.france-bioinformatique.fr/en/home/" class="text-white text-decoration-none" target="_blank" rel="noreferrer">
		French Institute of Bioinformatics
	</a>
</h5>

<a href="https://www.france-bioinformatique.fr/en/home/" target="_blank" rel="noreferrer">
	<img src="/logo-ifb.svg" alt="Logo of the IFB" width="250" />
</a>

<TutorialList
	categories={[
		{
			name: "",
			tutorials: [ifblinuxbasics1, ifblinuxbasics2, ifblinuxbasics3]
		}
	]}
/>

<p class="text-muted small">
	IFB tutorials are licensed under a
	<a href="https://creativecommons.org/licenses/by-sa/4.0/legalcode" target="_blank" rel="noreferrer"> CC-BY-SA </a>
	license.
</p>
