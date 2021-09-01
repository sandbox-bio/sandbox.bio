<script>
import { Icon, Tooltip } from "sveltestrap";
import { progress } from "./stores/config";

export let items = [];
export let title = "Tutorials";
export let max = Infinity;  // max number of tutorials to list
</script>

<div class="row mt-5">
	<h3 class="pb-2">{title}</h3>
</div>
<div class="row align-items-md-stretch">
	{#each items.slice(0, max) as info}
		<div class="col-md-6 col-lg-4 col-xl-3 mt-2">
			<div class="h-100 p-3 border rounded-3">
				{#each (info.difficulty || []) as tag}
					<span class="badge" class:bg-success={tag == "beginner"} class:bg-danger={tag == "difficult"} style={tag == "intermediate" ? "background-color:#fd7e14" : ""}>{tag}</span>
				{/each}
				{#each (info.tags || []) as tag}
					<span class="badge bg-primary me-1 mb-2">{tag}</span>
				{/each}
				<h4>
					{info.name}
					{#if info.id in $progress}
						{#if $progress[info.id].step == (info.steps.length - 1)}
							<Icon id={`icon-${info.id}`} class="float-end text-primary" name="check-circle" />
							<Tooltip target={`icon-${info.id}`}>
								Done!
							</Tooltip>
						{:else if $progress[info.id].step && $progress[info.id].step > 0}
							<a href={`/tutorials?id=${info.id}&step=${$progress[info.id].step}`}>
								<Icon id={`icon-${info.id}`} class="float-end text-primary" name="hourglass-split" />
							</a>	
							<Tooltip target={`icon-${info.id}`}>
								You're at step {$progress[info.id].step + 1} / {info.steps.length}
							</Tooltip>
						{/if}
					{/if}
				</h4>
				<p>{@html info.description}</p>
				<!-- "Explore" listings -->
				{#if info.url}
					<a href={info.url} class="btn btn-outline-primary text-center" target="_blank">Launch</a>
				<!-- "Tutorials" listings -->
				{:else}
					<a href={`/tutorials?id=${info.id}&step=${info.id in $progress ? $progress[info.id].step : 0}`} class="btn btn-outline-primary text-center">
						{#if info.id in $progress}
							Resume
						{:else}
							Launch
						{/if}
					</a>
				{/if}
			</div>
		</div>
	{/each}
</div>
