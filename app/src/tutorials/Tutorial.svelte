<script>
import { DropdownItem, Offcanvas } from "sveltestrap";
import { config, progress } from "./stores/config";
import { status } from "./stores/status";
import { tutorial, tutorials } from "./stores/tutorials";
import Terminal from "./terminal/Terminal.svelte";
import IDE from "./components/IDE.svelte";

export let id;
export let step = 0;

// State
$tutorial = $tutorials.find(t => t.id == id);
const tocToggle = () => tocOpen = !tocOpen;
let tocOpen = false;
let stepInfo = {};

// Reactive statements
$: nextStep(step);
$: nbSteps = $tutorial.steps.length;

function nextStep(step)
{
	stepInfo = $tutorial.steps[step];

	// Update URL
	const url = new URL(window.location);
	if(+url.searchParams.get("step") != +step) {
		url.searchParams.set("step", step);

		// Handle analytics before updating the URL
		try {
			fetch(`${$config.api}/ping`, {
				method: "POST",
				mode: "no-cors",
				body: JSON.stringify({
					tutorial: $tutorial.id,
					from: +new URL(window.location).searchParams.get("step"),
					to: +new URL(url).searchParams.get("step")
				})
			});
		} catch (error) {}

		window.history.pushState({}, "", url);
		$tutorial.step = step;
	}

	// Update progress in one shot (each time change $progress, makes call to DB)
	let progressNew = $progress || {};
	if(!($tutorial.id in progressNew))
		progressNew[$tutorial.id] = { step: 0 };
	// But only if the current step is greater!
	if(step > progressNew[$tutorial.id].step) {
		progressNew[$tutorial.id].step = step;
		$progress = progressNew;
	}

	// Scroll to the top when navigate pages
	if(document.getElementById("tutorial-sidebar"))
		document.getElementById("tutorial-sidebar").scrollTop = 0;
}

$tutorial.step = step;
</script>

<div class="container-fluid pb-3">
	<div class="d-grid gap-3" style="grid-template-columns: {nbSteps > 0 ? "1fr 2fr" : ""}; height:85vh; max-height:85vh">
		{#if $tutorial.steps.length > 0}
			<div class="bg-light border rounded-3 p-2 d-flex align-items-end flex-column">
				<div id="tutorial-sidebar" class="w-100 p-2 mb-auto" style="max-height:77vh; overflow-y:scroll; overflow-x:hidden">
					<h4>{stepInfo.name || $tutorial.name}</h4>
					{#if stepInfo.subtitle}
						<h6>{@html stepInfo.subtitle}</h6>
					{/if}
					{#if step == 0 && $tutorial.tags.length > 0}
						<div class="row mb-2">
							<h6>
								{#each $tutorial.tags as tag}
									<span class="badge bg-primary ms-1">
										{tag}
									</span>
								{/each}
								{#each $tutorial.difficulty as tag}
									<span class="badge" class:bg-success={tag == "beginner"} class:bg-danger={tag == "difficult"} style={tag == "intermediate" ? "background-color:#fd7e14" : ""}>{tag}</span>
								{/each}
							</h6>
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
						<div class="d-flex justify-content-between">
							<div>
								<button type="button" class="btn btn-sm" on:click={() => step--} class:btn-primary={step != 0} class:btn-secondary={step == 0} disabled={step == 0}>&larr; Previous</button>
								<button class="btn btn-sm" on:click={() => step++} class:btn-primary={step != $tutorial.steps.length - 1} class:btn-secondary={step == $tutorial.steps.length - 1} disabled={step == $tutorial.steps.length - 1}>Next &rarr;</button>	
							</div>
							<div>
								<a href="https://github.com/sandbox-bio/feedback/discussions/categories/questions" target="_blank">
									<span class="badge rounded-pill bg-secondary">Help</span>
								</a>
								<span on:click={tocToggle} class="badge rounded-pill bg-info">{step + 1} / {$tutorial.steps.length}</span>
							</div>
						</div>
					</div>
				</div>
			</div>
		{/if}
		{#if $tutorial.terminal !== false}
			<div id="terminal-wrapper" class="border rounded-3 p-2">
				<Terminal on:status={event => $status.terminal = event.detail} files={$tutorial.files} init={$tutorial.init} tools={$tutorial.tools} intro={$tutorial.intro} pwd={$tutorial.pwd} />
			</div>
		{/if}
		{#if $tutorial.ide === true}
			<IDE
				fn="dna_to_rna"
				code={`def dna_to_rna(t):\n\t# Your solution goes here\n\tpass\n`}
				input="GATGGAACTTGACTACGTAAATT"
				expectedInput="GATGGAACTTGACTACGTAAATT"
				expectedOutput="GAUGGAACUUGACUACGUAAAUU"
				/>
		{/if}
	</div>
</div>

<!-- TODO: this throws "Uncaught TypeError: $context is undefined" when click a lesson, but still seems to work -->
<Offcanvas header="Lessons" placement="end" isOpen={tocOpen} toggle={tocToggle}>
	<DropdownItem header>Introduction</DropdownItem>
	{#each $tutorial.steps as s, i}
		{#if s.header}
			<DropdownItem header><br />{s.name}</DropdownItem>
		{/if}
		<DropdownItem on:click={() => step = i}>
			{#if i == step}
				&rarr; <strong>{@html s.subtitle || s.name}</strong>
			{:else}
				<span style="visibility:hidden">&rarr;</span> {@html s.subtitle || s.name}
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
