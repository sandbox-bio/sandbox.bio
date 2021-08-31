<script>
import "bootstrap/dist/js/bootstrap.bundle";
import "bootstrap/dist/css/bootstrap.min.css";

import Home from "./routes/Home.svelte";
import Tutorial from "./tutorials/Tutorial.svelte";
import Terminal from "./terminal/Terminal.svelte";
import Listings from "./components/Listings.svelte";
import { config } from "./stores/config";
import { tutorials } from "./stores/tutorials";

// Config
const intro = $config.playground;
const path = window.location.pathname;
const params = new URL(window.location).searchParams;
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
			<div class="input-group mb-3 nav-pills">
				<a role="button" class="btn btn-link text-decoration-none pe-1 nav-link" class:active={path == "/tutorials"} href="/tutorials">Tutorials</a>
				<button type="button" class="btn btn-link text-decoration-none dropdown-toggle dropdown-toggle-split ps-2 pe-3 nav-link" class:active={path == "/tutorials"} data-bs-toggle="dropdown" aria-expanded="false">
					<span class="visually-hidden">Toggle Dropdown</span>
				</button>
				<ul class="dropdown-menu">
					{#each $tutorials as tutorial}
						<li><a class="dropdown-item" href="/tutorials?id={tutorial.id}">{tutorial.name}</a></li>
					{/each}
				</ul>
			</div>
		</li>
		<li class="nav-item">
			<a href="/playground" class="nav-link" class:active={path == "/playground"}>Playground</a>
		</li>
	</ul>
</header>

<main role="main" class="ps-4 pe-4">
	{#if path == "/"}
		<Home />
	{:else if path.startsWith("/tutorials")}
		{#if params.get("id")}
			<Tutorial id={params.get("id")} step={+params.get("step") || 0} />
		{:else}
			<Listings items={$tutorials} />
		{/if}
	{:else if path.startsWith("/playground")}
		<div class="p-2" style="background-color:#000">
			<Terminal {intro} files={$tutorials[1].files} />
		</div>
	{/if}
</main>

