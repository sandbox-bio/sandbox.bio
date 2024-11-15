<script>
import TutorialList from "$src/components/TutorialList.svelte";
import { config as datepro1 } from "$content/datepro-01-unix-shell/config";
import { config as datepro2 } from "$content/datepro-02-data-retrieval/config";
import { config as datepro3 } from "$content/datepro-03-data-extraction/config";
import { config as datepro4 } from "$content/datepro-04-task-repetition/config";
import { config as datepro5 } from "$content/datepro-05-xml-processing/config";
import { config as datepro6 } from "$content/datepro-06-text-retrieval/config";
import { config as datepro7 } from "$content/datepro-07-text-processing/config";
import { config as datepro8 } from "$content/datepro-08-semantic-processing/config";

let title = "Data and Text Processing";
let description = "Command-line tutorials on data and text processing, by Francisco M. Couto.";
let url = "https://sandbox.bio/community/data-processing";
let image = "https://sandbox.bio/assets/community/data-processing/screenshot.png";
</script>

<svelte:head>
	<title>Data and Text Processing for Health and Life Sciences - Community Tutorials - sandbox.bio</title>

	<meta name="og:type" content="website" />
	<meta name="og:title" content={title} />
	<meta name="og:description" content={description} />
	<meta name="og:url" content={url} />
	<meta name="og:image" content={image} />

	<meta name="twitter:title" content={title} />
	<meta name="twitter:description" content={description} />
	<meta name="twitter:url" content={url} />
	<meta name="twitter:image" content={image} />
</svelte:head>

<p class="m-0 text-muted opacity-75 fw-light small">COMMUNITY TUTORIALS</p>
<h4 class="my-0">Data and Text Processing for Health and Life Sciences</h4>
<p class="text-muted">
	By
	<a href="https://webpages.ciencias.ulisboa.pt/~fjcouto/" class="text-decoration-none" target="_blank" rel="noreferrer">Francisco M. Couto</a>
</p>

<a href="http://labs.rd.ciencias.ulisboa.pt/book/" target="_blank" rel="noreferrer">
	<img src="/assets/community/data-processing/logo.jpg" alt="Book cover" height="200" />
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
