<script>
import { Badge, Button, Card, Col, Icon, Input, Row, Tooltip } from "sveltestrap";
import { user } from "$stores/user";
import { progress } from "$stores/progress";
import MailingList from "./MailingList.svelte";

const TAG_COLORS = {
	beginner: "primary",
	difficult: "warning",
	advanced: "danger",
	intermediate: "warning"
};

export let categories = [];
</script>

{#each categories as category, i}
	<h5 class:mt-4={i > 0}>
		{#if category.icon}
			<Icon name={category.icon} class="text-primary" />
		{/if}
		{category.name}
	</h5>

	<Row cols={{ lg: 3, md: 2, sm: 1, xs: 1 }}>
		{#each category.tutorials as tutorial}
			{@const haveProgressInfo = $user?.email && tutorial.id in $progress}
			{@const numSteps = tutorial.steps?.length || 0}
			{@const currStep = (haveProgressInfo && Math.min($progress[tutorial.id].step, numSteps)) || -1}
			{@const isDone = currStep == numSteps - 1}
			{@const isInProgress = currStep > 0}
			{@const elementId = `tutorial-${tutorial.id}-${Date.now()}`}
			<!-- Otherwise when navigate from other pages, tooltip doesn't show up -->
			<Col class="my-2">
				<!-- Show progress tooltip, only for tutorials -->
				<Tooltip target={elementId}>
					{#if isDone}
						You completed this tutorial!
					{:else if isInProgress}
						You completed {currStep} / {numSteps} steps. Click to continue.
					{:else}
						Click to start this tutorial!
					{/if}
				</Tooltip>

				<!-- Tutorial card -->
				<Card id={elementId} class="h-100 listing-card p-3 bg-opacity-25 {isDone ? 'bg-success' : isInProgress ? 'bg-primary' : ''}">
					<!-- Tags -->
					<div>
						{#each tutorial.difficulty || [] as tag}
							<Badge color="{TAG_COLORS[tag]} bg-opacity-75">{tag}</Badge>
						{/each}
						{#each tutorial.tags || [] as tag}
							<!-- Show tags as primary if not listing tutorials (if tutorials, too many colors) -->
							{@const color = tutorial.url ? "primary" : "secondary"}
							<Badge color="{color} bg-opacity-75" class="me-1 mb-1">{tag}</Badge>
						{/each}
						{#if tutorial.new}
							<Badge color="danger bg-opacity-100">new</Badge>
						{/if}

						<!-- Icon to represent tutorial done or not? -->
						<span class="float-end">
							{#if isDone}
								<Icon name="check-circle-fill" class="text-success h4" />
							{:else if isInProgress}
								<Icon name="hourglass-split" class="text-primary h4" />
							{/if}
						</span>
					</div>

					<!-- Tutorial name and description -->
					<h5 class="mt-1">{tutorial.name}</h5>
					<p>{@html tutorial.description}</p>

					<!-- Launch link -->
					<div>
						{#if isDone || !isInProgress}
							<Button outline color={isDone ? "success" : "primary"} class="stretched-link" href="/tutorials/{tutorial.id}">Start</Button>
						{:else}
							<Button color="primary" class="stretched-link" href="/tutorials/{tutorial.id}/{$progress[tutorial.id].step}">Resume</Button>
						{/if}
					</div>
				</Card>
			</Col>
		{/each}
	</Row>

	{#if category.mailinglist}
		<MailingList />
	{/if}
{/each}

<style>
:global(.listing-card:hover) {
	background-color: #eee !important;
}
</style>
