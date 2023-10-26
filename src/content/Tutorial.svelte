<script>
import { onDestroy, onMount } from "svelte";
import { Button, DropdownItem, Offcanvas } from "sveltestrap";
import { afterNavigate } from "$app/navigation";
import { progress } from "$stores/config";
import { status } from "$stores/status";
import { tutorials } from "$stores/tutorials";
import { tutorial } from "$stores/tutorial";
import Alert from "$src/components/Alert.svelte";
import Terminal from "$components/Terminal.svelte";
import IGV from "$components/igv/IGV.svelte";
import IDE from "$components/ExerciseRosalind.svelte";

export let id;
export let step = 0;

// State
$tutorial = $tutorials.find((t) => t.id == id) || { steps: [] };
const tocToggle = () => (tocOpen = !tocOpen);
let tocOpen = false;
let stepInfo = {};
let rosalind = {};

// Reactive statements
$: $tutorial.step = step;
$: nextStep(step);
$: nbSteps = $tutorial.steps.length;
$: isFirstStep = step === 0;
$: isLastStep = step === $tutorial.steps.length - 1;
$: if ($tutorial.ide === true) {
	rosalind = $tutorial.steps[step].rosalind;
	rosalind.fn = `solution_${rosalind.id.toLowerCase()}`;
}

// Make sure to reset scroll position when navigating between steps
afterNavigate(() => {
	document.getElementById("tutorial-sidebar")?.scrollTo(0, 0);
});

// Handle analytics
async function logStep(from, to) {
	try {
		fetch(`/api/v1/ping`, {
			method: "POST",
			mode: "no-cors",
			body: JSON.stringify({
				from,
				to,
				tutorial: $tutorial.id
			})
		});
	} catch (error) {
		console.error(error);
	}
}

function nextStep(step) {
	if (!$tutorial.id) return;
	stepInfo = $tutorial.steps[step];

	// Update progress in one shot (each time change $progress, makes call to DB)
	let progressNew = $progress || {};
	if (!($tutorial.id in progressNew)) progressNew[$tutorial.id] = { step: 0 };
	// But only if the current step is greater!
	if (step > progressNew[$tutorial.id].step) {
		progressNew[$tutorial.id].step = step;
		$progress = progressNew;
	}
}

// When first visit a tutorial
onMount(() => {
	logStep(null, step);
});

// When leave tutorial, reset state
onDestroy(() => {
	$tutorial = {};
});
</script>

<div class="container-fluid pb-3 px-0">
	<div class="d-grid gap-2" style="grid-template-columns: {nbSteps > 0 ? '1fr 2fr' : ''}; height:85vh; max-height:85vh">
		{#if $tutorial.steps.length > 0}
			<div class="bg-light border rounded-3 p-2 d-flex align-items-end flex-column">
				<div id="tutorial-sidebar" class="w-100 p-2 mb-auto" style="max-height:77vh; overflow-y:scroll; overflow-x:hidden">
					<h4>{stepInfo.name || $tutorial.name}</h4>
					{#if stepInfo.subtitle || ($tutorial.step == 0 && $tutorial.subtitle)}
						<h6>{@html stepInfo.subtitle || $tutorial.subtitle}</h6>
					{/if}
					{#if step == 0 && $tutorial.tags.length > 0}
						<div class="row mb-2">
							<h6>
								{#each $tutorial.tags as tag}
									<span class="badge bg-primary ms-1" class:bg-info={tag == "course"}>
										{tag}
									</span>
								{/each}
								{#each $tutorial.difficulty as tag}
									<span
										class="badge"
										class:bg-success={tag == "beginner"}
										class:bg-danger={tag == "difficult"}
										style={tag == "intermediate" ? "background-color:#fd7e14" : ""}>{tag}</span
									>
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
								<Button
									href={isFirstStep ? "#" : `/tutorials/${$tutorial.id}` + (step === 1 ? "" : `/${step - 1}`)}
									color={isFirstStep ? "secondary" : "primary"}
									size="sm"
									on:click={() => logStep(step, step - 1)}
								>
									&larr;<span class="mobile-hide">&nbsp;Previous</span>
								</Button>

								<Button
									href={isLastStep ? "#" : `/tutorials/${$tutorial.id}/${step + 1}`}
									color={isLastStep ? "secondary" : "primary"}
									size="sm"
									on:click={() => logStep(step, step + 1)}
								>
									<span class="mobile-hide">Next&nbsp;</span>&rarr;
								</Button>
							</div>
							<div>
								<a href="https://github.com/sandbox-bio/sandbox.bio/discussions" target="_blank" rel="noreferrer">
									<span class="badge rounded-pill bg-secondary">Help</span>
								</a>
								<span on:click={tocToggle} class="badge rounded-pill bg-info">{step + 1} / {$tutorial.steps.length}</span>
							</div>
						</div>
					</div>
				</div>
			</div>
		{/if}
		{#if $tutorial.ide === true}
			<IDE
				fn={rosalind.fn}
				fnParams={rosalind.params}
				code={`def ${rosalind.fn}(${rosalind.params.join(", ")}):\n    # Your solution goes here\n    return "answer"\n`}
				input={rosalind.sample_data}
				expectedInput={rosalind.sample_data}
				expectedOutput={rosalind.sample_output}
			/>
		{:else if $tutorial.igv === true}
			{@const config = { ...$tutorial.igvConfig.default, ...$tutorial.igvConfig[step] }}
			<IGV options={config} />
		{:else if $tutorial.id != null}
			<div id="terminal-wrapper" class="border rounded-3 p-2">
				<Terminal
					on:status={(event) => ($status.terminal = event.detail)}
					assets={$tutorial.assets}
					files={$tutorial.files}
					tools={$tutorial.tools}
					intro={$tutorial.intro}
					init={$tutorial.init}
				/>
			</div>
		{:else}
			<div>
				<Alert color="warning">
					<p>
						<strong>This tutorial does not exist.</strong>
					</p>
					<p>
						Please <a href="https://github.com/sandbox-bio/sandbox.bio/discussions/new?category=general">reach out to us</a> and let us know how you got
						here so we can fix the broken link. Please include the original page that brought you here, thank you!
					</p>
				</Alert>
				<Button color="primary" href="/tutorials">Browse tutorials &rarr;</Button>
			</div>
		{/if}
	</div>
</div>

<!-- Note: this throws "Uncaught TypeError: Cannot read properties of undefined (reading 'autoClose')"
	when clicking a lesson, but it still works -->
<Offcanvas header="Lessons" placement="end" isOpen={tocOpen} toggle={tocToggle}>
	<DropdownItem header class="text-muted">Introduction</DropdownItem>

	{#each $tutorial.steps as s, i}
		{#if s.header}
			<DropdownItem header class="mt-4 text-muted">
				{s.name}
			</DropdownItem>
		{/if}

		<DropdownItem href="/tutorials/{$tutorial.id}/{i}" class="my-2">
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

@media only screen and (max-width: 768px) {
	.mobile-hide {
		display: none;
	}
}
</style>
