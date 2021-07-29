<script>
import Home from "./Home.svelte";
import Tutorials from "./Tutorials.svelte";
import "bootstrap/dist/js/bootstrap.bundle";
import "bootstrap/dist/css/bootstrap.min.css";

// State
const path = window.location.pathname;
</script>

<header class="d-flex flex-wrap justify-content-center py-3 border-bottom ms-4 me-4 mb-3">
	<a href="/" class="d-flex align-items-center mb-md-0 me-md-auto text-dark text-decoration-none">
		<span class="fs-4">🧬 sandbox.bio</span>
	</a>

	<ul class="nav nav-pills">
		<li class="nav-item">
			<a href="/playground" class="nav-link" class:active={path == "/playground"}>Playground</a>
		</li>
		<li class="nav-item dropdown">
			<!-- svelte-ignore a11y-no-onchange -->
			<a class="nav-link dropdown-toggle" class:active={path == "/playground"} href="#" id="navTutorials" role="button" data-bs-toggle="dropdown" aria-expanded="false">
				Tutorials
			</a>
			<ul class="dropdown-menu" aria-labelledby="navTutorials">
				<li><a class="dropdown-item" href="/tutorials">Introduction to bedtools</a></li>
				<li><a class="dropdown-item" href="/tutorials">Introduction to bowtie2</a></li>
			</ul>
		</li>
		<li class="nav-item dropdown">
			<!-- svelte-ignore a11y-no-onchange -->
			<a class="nav-link dropdown-toggle" href="#" id="navExplore" role="button" data-bs-toggle="dropdown" aria-expanded="false">
				Explore
			</a>
			<ul class="dropdown-menu" aria-labelledby="navExplore">
				<li><a class="dropdown-item" href="https://alignment.sandbox.bio" target="_blank">Sequence alignment</a></li>
				<li><a class="dropdown-item" href="https://tsne.sandbox.bio" target="_blank">tSNE algorithm</a></li>
				<li><a class="dropdown-item" href="https://fastq.sandbox.bio" target="_blank">FASTQ QC metrics</a></li>
			</ul>
		</li>
	</ul>
</header>

<main role="main" class="ps-4 pe-4">
	{#if path == "/"}
		<Home />
	{:else if path == "/tutorials"}
		<Tutorials />
	{/if}
</main>
