<script>
import TutorialList from "$src/components/TutorialList.svelte";
import { config as ifblinuxbasics1 } from "$content/ifb-linux-basics-1/config.js";
import { config as ifblinuxbasics2 } from "$content/ifb-linux-basics-2/config.js";
import { config as ifblinuxbasics3 } from "$content/ifb-linux-basics-3/config.js";

import { config as datepro1 } from "$content/datepro-01-unix-shell/config";
import { config as datepro2 } from "$content/datepro-02-data-retrieval/config";
import { config as datepro3 } from "$content/datepro-03-data-extraction/config";
import { config as datepro4 } from "$content/datepro-04-task-repetition/config";
import { config as datepro5 } from "$content/datepro-05-xml-processing/config";
import { config as datepro6 } from "$content/datepro-06-text-retrieval/config";
import { config as datepro7 } from "$content/datepro-07-text-processing/config";
import { config as datepro8 } from "$content/datepro-08-semantic-processing/config";
</script>

<svelte:head>
	<title>Data and Text Processing for Health and Life Sciences - sandbox.bio</title>
</svelte:head>

<p class="m-0 text-muted opacity-75 fw-light small">COMMUNITY TUTORIALS</p>
<h4 class="my-0">Data and Text Processing for Health and Life Sciences</h4>
<p class="text-muted">
	By
	<a href="https://webpages.ciencias.ulisboa.pt/~fjcouto/" class="text-decoration-none" target="_blank" rel="noreferrer"> Francisco M. Couto </a>
</p>

<a href="http://labs.rd.ciencias.ulisboa.pt/book/" target="_blank" rel="noreferrer">
	<img src="/logo-book-data-processing.jpg" alt="Book cover" height="200" />
</a>

<TutorialList
	categories={[
		{
			name: "",
			tutorials: [datepro1, datepro2, datepro3, datepro4, datepro5, datepro6, datepro7, datepro8]
		}
	]}
/>

<p class="text-muted small">
	<i>Data and Text Processing for Health and Life Sciences</i> tutorials are licensed under a
	<a href="https://creativecommons.org/licenses/by/4.0/" target="_blank" rel="noreferrer">CC BY 4.0</a>
	license.
</p>
