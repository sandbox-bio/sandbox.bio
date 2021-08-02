<script>
import "bootstrap/dist/js/bootstrap.bundle";
import "bootstrap/dist/css/bootstrap.min.css";

import Home from "routes/Home.svelte";
import Tutorials from "routes/Tutorials.svelte";
import Terminal from "./Terminal.svelte";
import { config } from "config";

// State
const path = window.location.pathname;
const intro = `# This playground is for open-ended exploration.
# For guided tutorials, see https://sandbox.bio/tutorials.
#
# Example:
#   samtools view -o test.bam /samtools/examples/toy.sam
#   samtools index test.bam
#   ls test.bam.bai
#   samtools idxstats test.bam
`;
</script>

<svelte:head>
	<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.4.1/font/bootstrap-icons.css">
</svelte:head>

<header class="d-flex flex-wrap justify-content-center py-3 border-bottom ms-4 me-4 mb-3">
	<a href="/" class="d-flex align-items-center mb-md-0 me-md-auto text-dark text-decoration-none">
		<span class="fs-4">🧬 sandbox.bio</span>
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
	</ul>
</header>

<main role="main" class="ps-4 pe-4">
	{#if path == "/"}
		<Home />
	{:else if path == "/tutorials"}
		<Tutorials />
	{:else if path == "/playground"}
		<div class="p-2" style="background-color:#000">
			<Terminal {intro} files={$config.tutorials[0].files} />
		</div>
	{/if}
</main>
