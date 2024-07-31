<script>
import { Badge, Button, Card, Col, Row } from "sveltestrap";
export let items = [];

const tagColors = {
	beginner: "primary",
	difficult: "danger",
	intermediate: "warning"
};
</script>

<Row cols={{ lg: 3, md: 2, sm: 1, xs: 1 }}>
	{#each items as info}
		<Col class="my-2">
			<Card class="h-100 listing-card p-3">
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
				</div>

				<h4>{info.name}</h4>
				<p>{@html info.description}</p>

				<div><Button outline color="primary" class="stretched-link" href={info.url}>Launch</Button></div>
			</Card>
		</Col>
	{/each}
</Row>
