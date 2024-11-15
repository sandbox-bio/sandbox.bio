<script>
import TutorialList from "$src/components/TutorialList.svelte";
import { config as ifblinuxbasics1 } from "$content/ifb-linux-basics-1/config.js";
import { config as ifblinuxbasics2 } from "$content/ifb-linux-basics-2/config.js";
import { config as ifblinuxbasics3 } from "$content/ifb-linux-basics-3/config.js";
</script>

<p class="m-0 text-muted opacity-75 fw-light small">COMMUNITY TUTORIALS</p>
<h4 class="my-0">Linux Basics</h4>
<p class="text-muted">
	By the
	<a href="https://www.france-bioinformatique.fr/en/home/" class="text-decoration-none" target="_blank" rel="noreferrer">
		French Institute of Bioinformatics
	</a>
</p>

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
