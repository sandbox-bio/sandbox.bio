<script>
import { Badge, Button, Icon } from "sveltestrap";
import { user } from "$stores/user";
import { progress } from "$stores/progress";

export let items = [];
export let title = "Tutorials";
export let colMd = 6; // == 12 / how many boxes we can fit on M screens
export let colLg = 4; // == 12 / how many boxes we can fit on L screens
export let colXxl = 4; // == 12 / how many boxes we can fit on XXL screens
export let skip = []; // Specific tutorial IDs to not show
export let showUnlisted = false;

const tagColors = {
	beginner: "primary",
	difficult: "danger",
	intermediate: "warning"
};
</script>

{#if title}
	<h3>{title}</h3>
{/if}

<div class="row align-items-md-stretch">
	{#each items.filter((t) => !skip.includes(t.id) && (showUnlisted || t.listed !== false) && (t.steps?.length > 0 || t.url)) as info, i}
		{@const haveProgressInfo = $user && info.id in $progress}
		{@const isDone = haveProgressInfo && $progress[info.id].step == info.steps.length - 1}
		{@const isInProgress = haveProgressInfo && $progress[info.id].step && $progress[info.id].step > 0}

		{#if info.divider}
			<h5 class:mt-4={i > 0}>{info.divider}</h5>
		{/if}

		<div id="tutorial-{info.id}" class="col-md-{colMd} col-lg-{colLg} col-xxl-{colXxl} mt-2">
			<div
				class="listing-card h-100 p-3 border rounded-3 position-relative d-flex flex-column"
				class:bg-success={isDone}
				class:bg-opacity-25={isDone}
			>
				<!-- Tags -->
				<div>
					{#each info.difficulty || [] as tag}
						<Badge color="{tagColors[tag]} bg-opacity-75">{tag}</Badge>
					{/each}
					{#each info.tags || [] as tag}
						<!-- Show tags as primary if not listing tutorials (if tutorials, too many colors) -->
						{@const color = info.url ? "primary" : "secondary"}
						<Badge color="{color} bg-opacity-75" class="me-1 mb-1">{tag}</Badge>
					{/each}

					<span class="float-end">
						{#if isDone}
							<Icon name="check-circle-fill" class="text-success h4" />
						{:else if isInProgress}
							<Icon name="hourglass-split" class="text-primary h4" />
						{/if}
					</span>
				</div>

				<!-- Tutorial Info -->
				<h4>{info.name}</h4>
				<p>{@html info.description}</p>

				<!-- Launch link -->
				<div>
					<!-- "Explore" listings -->
					{#if info.url}
						<Button color="primary" class="stretched-link" href={info.url}>Launch</Button>

						<!-- "Tutorials" listings -->
					{:else if haveProgressInfo}
						{#if isDone}
							<Button color="success" class="stretched-link" href="/tutorials/{info.id}" outline>Launch</Button>
						{:else}
							<Button color="primary" class="stretched-link" href="/tutorials/{info.id}/{$progress[info.id].step}">Resume</Button>
						{/if}
					{:else}
						<Button color="primary" class="stretched-link" href="/tutorials/{info.id}">Launch</Button>
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
