<script>
import { Button, Icon } from "sveltestrap";
import { categories, explore, playgrounds } from "$stores/tutorials";
import Listings from "$components/Listings.svelte";
import TutorialList from "$src/components/TutorialList.svelte";
</script>

<svelte:head>
	<title>sandbox.bio</title>
</svelte:head>

<div class="container-fluid pt-3">
	<div class="row pb-0 pe-lg-0 bg-light align-items-center rounded-3 border shadow-lg">
		<div class="col-lg-6 p-3 p-lg-4 pt-lg-3">
			<h1 class="fw-bold" style="font-size:35px">Interactive bioinformatics tutorials</h1>
			<p class="lead my-4">Learn bioinformatics from your browser.<br />Everything runs in a sandbox, so you can experiment all you want.</p>
			<Button size="lg" color="primary" href="/tutorials">Get started <Icon name="arrow-right" /></Button>
			<Button size="lg" color="outline-primary" href="/about">Learn more</Button>
		</div>
		<div class="align-center col-lg-5 p-0 offset-lg-1 overflow-hidden shadow-lg">
			<img class="rounded-lg-3" src="/screenshot-cli.png" alt="Screenshot of terminal and exercises" width="800" />
		</div>
	</div>
</div>

<h4 class="mt-5">Tutorials</h4>
<TutorialList categories={$categories} />

<h4 class="mt-5">Playgrounds</h4>
<Listings items={$playgrounds} />

<h4 class="mt-5">Explore</h4>
<Listings items={$explore} />

<style>
.bg-light {
	background-color: rgb(232, 240, 240) !important;
}
</style>
