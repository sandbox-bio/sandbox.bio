<script>
import TutorialList from "$src/components/TutorialList.svelte";
import { config as ifblinuxbasics1 } from "$content/ifb-linux-basics-1/config.js";
import { config as ifblinuxbasics2 } from "$content/ifb-linux-basics-2/config.js";
import { config as ifblinuxbasics3 } from "$content/ifb-linux-basics-3/config.js";

import { config as datepro1 } from "$content/datepro-01-unix-shell/config";
import { config as datepro2 } from "$content/datepro-02-data-retrieval/config";
import { config as datepro3 } from "$content/datepro-03-data-extraction/config";
</script>

<svelte:head>
	<title>Community - sandbox.bio</title>
</svelte:head>

<h3 class="my-3">Community Tutorials</h3>

<h5 class="bg-primary-subtle p-2">
	<a href="https://www.france-bioinformatique.fr/en/home/" class="text-decoration-none" target="_blank" rel="noreferrer">
		French Institute of Bioinformatics
	</a>
</h5>

<a href="https://www.france-bioinformatique.fr/en/home/" target="_blank" rel="noreferrer">
	<img src="/logo-ifb.svg" alt="IFB logo" height="100" />
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
	<i>IFB</i> tutorials are licensed under a
	<a href="https://creativecommons.org/licenses/by-sa/4.0/legalcode" target="_blank" rel="noreferrer"> CC-BY-SA </a>
	license.
</p>

<h5 class="bg-primary-subtle p-2 mt-5">
	<a href="http://labs.rd.ciencias.ulisboa.pt/book/" class="text-decoration-none" target="_blank" rel="noreferrer">
		Data and Text Processing for Health and Life Sciences, by Francisco M. Couto
	</a>
</h5>

<a href="http://labs.rd.ciencias.ulisboa.pt/book/" target="_blank" rel="noreferrer">
	<img
		src="https://webpages.ciencias.ulisboa.pt/~fjcouto/Data%20and%20Text%20Processing%20for%20Health%20and%20Life%20Sciences.jpg"
		alt="Book cover"
		height="150"
	/>
</a>

<TutorialList
	categories={[
		{
			name: "",
			tutorials: [datepro1, datepro2, datepro3]
		}
	]}
/>

<p class="text-muted small">
	<i>Data and Text Processing for Health and Life Sciences</i> tutorials are licensed under a
	<a href="https://creativecommons.org/licenses/by/4.0/" target="_blank" rel="noreferrer">CC BY 4.0</a>
	license.
</p>
