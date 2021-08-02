<script>
import { config } from "config";
import Terminal from "./Terminal.svelte";
import { DropdownItem, Offcanvas } from "sveltestrap";

export let id;
export let step = 0;

// State
const tutorial = $config.tutorials.find(t => t.id == id);
const tocToggle = () => tocOpen = !tocOpen;
let tocOpen = false;
let stepInfo = {};

// Reactive statements
$: nextStep(step);

function nextStep(step)
{
	stepInfo = tutorial.steps[step];

	// Update URL
	const url = new URL(window.location);
	if(+url.searchParams.get("step") != +step) {
		url.searchParams.set("step", step);
		window.history.pushState({}, "", url);
	}

	// Scroll to the top when navigate pages
	if(document.getElementById("tutorial-sidebar"))
		document.getElementById("tutorial-sidebar").scrollTop = 0;
}
</script>

<div class="container-fluid pb-3">
	<div class="d-grid gap-3" style="grid-template-columns: 1fr 2fr; height:85vh; max-height:85vh">
		<div class="bg-light border rounded-3 p-2 d-flex align-items-end flex-column">
			<div id="tutorial-sidebar" class="w-100 p-2 mb-auto" style="max-height:77vh; overflow-y:scroll; overflow-x:hidden">
				<h4>{stepInfo.name || tutorial.name}</h4>
				{#if stepInfo.subtitle}
					<h6>{stepInfo.subtitle}</h6>
				{/if}
				{#if step == 0}
					<div class="row mb-2">
						<h6>
							{#each tutorial.tools as tag}
								<span class="badge bg-primary">{tag}</span>
							{/each}
							{#each tutorial.difficulty as tag}
								<span class="badge" class:bg-success={tag == "beginner"} class:bg-warning={tag == "intermediate"} class:bg-danger={tag == "difficult"}>{tag}</span>
							{/each}
						</h6>
						{#if tutorial.adapted_from}
							<span>Adapted from <a href={tutorial.adapted_from.link} target="_blank">{tutorial.adapted_from.name}</a></span>
						{/if}
					</div>
				{/if}
				<hr class="border-2 border-top border-secondary" />

				<div id="tutorial-wrapper" class="row" style="overflow-x: hidden">
					<div class="container">
						<svelte:component this={stepInfo.component} />
					</div>
				</div>
			</div>

			<div class="w-100 p-2 border-top pt-4">
				<div class="row">
					<div class="col-10">
						<button type="button" class="btn btn-sm" on:click={() => step--} class:btn-primary={step != 0} class:btn-secondary={step == 0} disabled={step == 0}>&larr; Previous</button>
						<button class="btn btn-sm" on:click={() => step++} class:btn-primary={step != tutorial.steps.length - 1} class:btn-secondary={step == tutorial.steps.length - 1} disabled={step == tutorial.steps.length - 1}>Next &rarr;</button>
					</div>
					<div class="col-2 text-end">
						<span on:click={tocToggle} class="badge rounded-pill bg-info">{step + 1} / {tutorial.steps.length}</span>
					</div>
				</div>
			</div>
		</div>
		<div id="terminal-wrapper" class="border rounded-3 p-2">
			<Terminal files={tutorial.files} />
		</div>
	</div>
</div>

<!-- TODO: this throws "Uncaught TypeError: $context is undefined" when click a lesson, but still seems to work -->
<Offcanvas header="Lessons" placement="end" isOpen={tocOpen} toggle={tocToggle}>
	<DropdownItem header>Introduction</DropdownItem>
	{#each tutorial.steps as s, i}
		{#if s.header}
			<DropdownItem header><br />{s.name}</DropdownItem>
		{/if}
		<DropdownItem on:click={() => step = i}>
			{#if i == step}
				&rarr; <strong>{s.subtitle || s.name}</strong>
			{:else}
				<span style="visibility:hidden">&rarr;</span> {s.subtitle || s.name}
			{/if}
		</DropdownItem>
	{/each}
</Offcanvas>

<style>
#terminal-wrapper {
	background-color: black;
}

.rounded-pill:hover {
	cursor: pointer;
}
</style>
