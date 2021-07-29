<script>
import { config } from "../tutorials/config";

export let maxListings = Infinity;  // max number of tutorials to list

let id = new URL(window.location).searchParams.get("id");
let step = new URL(window.location).searchParams.get("step");
</script>

<!-- List tutorials if none specified -->
{#if !id}
	<div class="row mt-5">
		<h3 class="pb-2">Tutorials</h3>
	</div>
	<div class="row align-items-md-stretch">
		{#each config.tutorials.slice(0, maxListings) as tutorial}
			<div class="col-md-3">
				<div class="h-100 p-5 border rounded-3">
					<h4>{tutorial.name}</h4>
					<p>{@html tutorial.description}</p>
					<a href="/tutorials?id={tutorial.id}" class="btn btn-outline-dark text-center">Launch</a>
				</div>
			</div>
		{/each}
	</div>
{/if}
