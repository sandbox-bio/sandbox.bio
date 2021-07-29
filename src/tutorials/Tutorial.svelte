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

				<hr />
			
			{:else}

				<h4>Step {step}</h4>
				

			{/if}


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
</style>
