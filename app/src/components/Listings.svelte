<script>
export let items = [];
export let title = "Tutorials";
export let max = Infinity;  // max number of tutorials to list
</script>

<div class="row mt-5">
	<h3 class="pb-2">{title}</h3>
</div>
<div class="row align-items-md-stretch">
	{#each items.slice(0, max) as info}
		<div class="col-md-6 col-lg-4 mt-2">
			<div class="h-100 p-3 border rounded-3">
				{#each (info.difficulty || []) as tag}
					<span class="badge" class:bg-success={tag == "beginner"} class:bg-danger={tag == "difficult"} style={tag == "intermediate" ? "background-color:#fd7e14" : ""}>{tag}</span>
				{/each}
				{#each (info.tags || []) as tag}
					<span class="badge bg-primary me-1 mb-2">{tag}</span>
				{/each}
				<h4>{info.name}</h4>
				<p>{@html info.description}</p>
				<a href={info.url || `/tutorials?id=${info.id}`} class="btn btn-outline-primary text-center">Launch</a>
			</div>
		</div>
	{/each}
</div>
