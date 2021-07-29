<script>
import Terminal from "../terminal/Terminal.svelte";
import { config } from "../config";

export let id;
export let step = 0;

const tutorial = config.tutorials.find(t => t.id == id);
</script>

<div class="container-fluid pb-3">
	<div class="d-grid gap-3" style="grid-template-columns: 1fr 3fr;">
		<div class="bg-light border rounded-3 p-2">
			{#if step == 0}
			<div class="row mb-2">
				<h4>{tutorial.name}</h4>
				<h6>
						{#each tutorial.tools as tag}
							<span class="badge bg-primary">{tag}</span>
						{/each}
						{#each tutorial.difficulty as tag}
							<span class="badge" class:bg-success={tag == "beginner"} class:bg-warning={tag == "intermediate"} class:bg-danger={tag == "difficult"}>{tag}</span>
						{/each}
					</h6>
					<h6>by <a href={tutorial.author.link} target="_blank">{tutorial.author.name}</a></h6>
				</div>
				<hr class="border-2 border-top border-secondary" />
			{/if}

			<div class="row" id="tutorial-wrapper">
				<svelte:component this={tutorial.steps[step].component} />
			</div>

			<hr class="border-2 border-top border-secondary" />
			<div class="row mt-4">
				<div class="col-10">
					<button type="button" class="btn btn-sm btn-secondary" on:click={() => step--} disabled={step == 0}>&larr; Previous</button>
					<button class="btn btn-sm btn-primary" on:click={() => step++} disabled={step == tutorial.steps.length - 1}>Next &rarr;</button>
				</div>
				<div class="col-2 text-end">
					<span class="badge rounded-pill bg-info">{step + 1} / {tutorial.steps.length}</span>
				</div>
			</div>
		</div>
		<div id="terminal-wrapper" class="border rounded-3 p-2">
			<Terminal />
		</div>
	</div>
</div>

<style>
#terminal-wrapper {
	background-color: black;
}

#tutorial-wrapper {
	font-size: 1.2em;
}
</style>
