<script>
import { Icon, Tooltip } from "sveltestrap";
import { progress } from "$stores/config";

export let items = [];
export let title = "Tutorials";
export let max = Infinity; // max number of tutorials to list
export let colMd = 6; // == 12 / how many boxes we can fit on medium screens
export let colLg = 4; // == 12 / how many boxes we can fit on large screens
export let colXxl = 3; // == 12 / how many boxes we can fit on xxl screens
export let skip = []; // Specific tutorial IDs to not show
</script>

{#if title}
	<div class="row mt-5">
		<h3 class="pb-2">{title}</h3>
	</div>
{/if}

<div class="row align-items-md-stretch">
	{#each items.slice(0, max).filter((t) => !skip.includes(t.id) && t.listed !== false && (t.steps?.length > 0 || t.url)) as info, i}
		{#if info.divider}
			<h5 class:mt-4={i > 0}>{info.divider}</h5>
		{/if}

		<div class="col-md-{colMd} col-lg-{colLg} col-xxl-{colXxl} mt-2">
			<div class="listing-card h-100 p-3 border rounded-3 position-relative d-flex flex-column">
				<!-- Tags -->
				<div>
					{#each info.difficulty || [] as tag}
						<span
							class="badge"
							class:bg-success={tag == "beginner"}
							class:bg-danger={tag == "difficult"}
							style={tag == "intermediate" ? "background-color:#fd7e14" : ""}>{tag}</span
						>
					{/each}
					{#each info.tags || [] as tag}
						<span class="badge bg-primary me-1 mb-2">{tag}</span>
					{/each}
				</div>

				<!-- Tutorial Info -->
				<h4>
					{info.name}
					{#if info.id in $progress}
						{#if $progress[info.id].step == info.steps.length - 1}
							<Icon id={`icon-${info.id}`} class="float-end text-success" name="check-circle-fill" />
							<Tooltip target={`icon-${info.id}`}>Done!</Tooltip>
						{:else if $progress[info.id].step && $progress[info.id].step > 0}
							<a href={`/tutorials/${info.id}/${$progress[info.id].step}`}>
								<Icon id={`icon-${info.id}`} class="float-end text-primary" name="circle-half" />
							</a>
							<Tooltip target={`icon-${info.id}`}>
								You're at step {$progress[info.id].step + 1} / {info.steps.length}
							</Tooltip>
						{/if}
					{/if}
				</h4>
				<p>{@html info.description}</p>

				<!-- Launch link -->
				<div>
					<!-- "Explore" listings -->
					{#if info.url}
						<a href={info.url} class="btn btn-outline-primary text-center mt-auto stretched-link">Launch</a>

						<!-- "Tutorials" listings -->
					{:else if info.id in $progress}
						{#if $progress[info.id].step == info.steps.length - 1}
							<a href={`/tutorials/${info.id}`} class="btn btn-outline-success text-center mt-auto stretched-link">Launch</a>
						{:else}
							<a href={`/tutorials/${info.id}/${$progress[info.id].step}`} class="btn btn-primary text-center mt-auto stretched-link"
								>Resume</a
							>
						{/if}
					{:else}
						<a href={`/tutorials/${info.id}`} class="btn btn-outline-primary text-center mt-auto stretched-link">Launch</a>
					{/if}
				</div>
			</div>
		</div>
	{/each}
</div>

<style>
.listing-card:hover {
	background-color: #eee !important;
}
</style>
