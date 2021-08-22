<script>
import { Modal } from "sveltestrap";
import "bootstrap/dist/js/bootstrap.bundle";
import "bootstrap/dist/css/bootstrap.min.css";

import Home from "./routes/Home.svelte";
import Tutorials from "./routes/Tutorials.svelte";
import Terminal from "./Terminal.svelte";
import { config } from "config";

// State
let aboutIsOpen = false;  // Whether "About" modal is showing or not
const path = window.location.pathname;
const intro = `\u001b[0;37m# This playground is for open-ended exploration.
# For guided tutorials, see https://sandbox.bio/tutorials
#
# Example:
#   samtools view -o test.bam /samtools/examples/toy.sam
#   samtools index test.bam
#   ls test.bam.bai
#   samtools idxstats test.bam  # idxstats uses the .bai file
\u001b[0m
`;
</script>

<svelte:head>
	<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.4.1/font/bootstrap-icons.css">
</svelte:head>

<header class="d-flex flex-wrap justify-content-center py-3 border-bottom ms-4 me-4 mb-3">
	<a href="/" class="d-flex align-items-center mb-md-0 me-md-auto text-dark text-decoration-none">
		<span class="fs-4">🧬 sandbox.bio</span>&nbsp;
	</a>

	<ul class="nav nav-pills">
		<li class="nav-item dropdown">
			<!-- svelte-ignore a11y-invalid-attribute -->
			<a href="#" class="nav-link dropdown-toggle" class:active={path == "/tutorials"} id="navTutorials" role="button" data-bs-toggle="dropdown" aria-expanded="false">
				Tutorials
			</a>
			<ul class="dropdown-menu" aria-labelledby="navTutorials">
				{#each $config.tutorials as tutorial}
					<li><a class="dropdown-item" href="/tutorials?id={tutorial.id}">{tutorial.name}</a></li>
				{/each}
			</ul>
		</li>
		<li class="nav-item dropdown">
			<!-- svelte-ignore a11y-invalid-attribute -->
			<a class="nav-link dropdown-toggle" href="#" id="navExplore" role="button" data-bs-toggle="dropdown" aria-expanded="false">
				Explore
			</a>
			<ul class="dropdown-menu" aria-labelledby="navExplore">
				<li><a class="dropdown-item" href="https://alignment.sandbox.bio" target="_blank">Sequence alignment</a></li>
				<li><a class="dropdown-item" href="https://tsne.sandbox.bio" target="_blank">tSNE algorithm</a></li>
				<li><a class="dropdown-item" href="https://fastq.sandbox.bio" target="_blank">FASTQ QC metrics</a></li>
				<li><a class="dropdown-item" href="https://wgsim.sandbox.bio" target="_blank">Simulate DNA sequences</a></li>
			</ul>
		</li>
		<li class="nav-item">
			<a href="/playground" class="nav-link" class:active={path == "/playground"}>Playground</a>
		</li>
		<li class="nav-item">
			<button class="btn btn-link text-decoration-none" style="padding-top:7px" on:click={() => aboutIsOpen = !aboutIsOpen}>About</button>
		</li>
	</ul>
</header>

<main role="main" class="ps-4 pe-4">
	{#if path == "/"}
		<Home />
	{:else if path.startsWith("/tutorials")}
		<Tutorials />
	{:else if path.startsWith("/playground")}
		<div class="p-2" style="background-color:#000">
			<Terminal {intro} files={$config.tutorials[0].files} />
		</div>
	{/if}
</main>

<Modal body scrollable size="lg" header="About sandbox.bio" toggle={() => aboutIsOpen = !aboutIsOpen} isOpen={aboutIsOpen}>
	<p>Interactive tutorials for exploring bioinformatics command-line tools in a secure sandbox.</p>

	<p class="lead fw-bold mt-4 mb-1">How does it work?</p>
	<p>We compiled commonly-used bioinformatics tools to WebAssembly so that they can run in your browser (see our <a href="https://github.com/biowasm/biowasm" target="_blank">biowasm</a> project).</p>

	<p>To simulate a terminal environment, we implemented features such as piping (<code>|</code>), file redirection (<code>></code>, <code>>></code>), process substitution (<code>&lt;()</code>), conditional commands (<code>&&</code>, <code>||</code>), asynchronous commands (<code>&</code>), variables (<code>abc=123</code>), and autocomplete. The terminal UI is based on <a href="https://github.com/xtermjs/xterm.js/" target="_blank">xterm.js</a>.</p>

	<p>GNU Coreutils such as <code>ls</code>, <code>cat</code>, <code>grep</code>, <code>head</code>, <code>wc</code>, and <code>echo</code> were also implemented, though keep in mind we don't currently support most flags (to address that in the future, we'll compile coreutils to WebAssembly!).
		All the files you read and write to are temporarily stored in memory using Emscripten's <a href="https://emscripten.org/docs/api_reference/Filesystem-API.html" target="_blank">virtual file system</a>.</p>

	<p class="lead fw-bold mt-4 mb-1">How to contribute</p>
	<p>If you have feedback, ideas for new tutorials, or if one of your own bioinformatics tutorials could benefit from being interactive, please <a href="mailto:robert.aboukhalil+sandboxbio@gmail.com">reach out</a>! Note that currently, only the C/C++ bioinformatics tools listed <a href="https://github.com/biowasm/biowasm#supported-tools" target="_blank">here</a> are supported.</p>

	<p class="lead fw-bold mt-4 mb-1">Author</p>
	<p>Built by <a href="https://www.robertaboukhalil.com/" target="_blank">Robert Aboukhalil</a>.</p>
</Modal>
