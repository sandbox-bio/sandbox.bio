<script>
import { config } from "config";
import Tutorial from "./Tutorial.svelte";

export let maxListings = Infinity;  // max number of tutorials to list

const id = new URL(window.location).searchParams.get("id");
const step = +new URL(window.location).searchParams.get("step") || 0;
const tutorial = $config.tutorials.find(t => t.id == id);
</script>

<!-- Load tutorial if we know which one we want -->
{#if tutorial}
	<Tutorial {id} {step} />

<!-- List tutorials if none specified -->
{:else}
	<div class="row mt-5">
		<h3 class="pb-2">Tutorials</h3>
	</div>
	<div class="row align-items-md-stretch">
		{#each $config.tutorials.slice(0, maxListings) as info}
			<div class="col-md-3">
				<div class="h-100 p-5 border rounded-3">
					<h4>{info.name}</h4>
					<p>{@html info.description}</p>
					<a href="/tutorials?id={info.id}" class="btn btn-outline-dark text-center">Launch</a>
				</div>
			</div>
		{/each}
	</div>
{/if}
