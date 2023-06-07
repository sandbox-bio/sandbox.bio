<script>
import { tutorials, explore, playgrounds } from "./stores/tutorials";
import Listings from "./components/Listings.svelte";
import Footer from "./components/Footer.svelte";
</script>

<div class="container-fluid pt-3">
	<div class="row pb-0 pe-lg-0 bg-light align-items-center rounded-3 border shadow-lg">
		<div class="col-lg-6 p-3 p-lg-4 pt-lg-3">
			<h1 class="display-6 fw-bold lh-1">Interactive bioinformatics tutorials</h1>
			<br />
			<p class="lead">Learn how to use bioinformatics tools right from your browser.<br />Everything runs in a sandbox, so you can experiment all you want.</p>
			<div class="d-grid gap-2 d-md-flex justify-content-md-start mb-4 mt-5 mb-lg-3">
				<a href="/tutorials?id=terminal-basics" class="btn btn-primary btn-lg px-4 me-md-2 fw-bold">Get started &rarr;</a>
			</div>
		</div>
		<div class="align-center col-lg-5 p-0 offset-lg-1 overflow-hidden shadow-lg">
			<img class="rounded-lg-3" src="images/cli.png" alt="Screenshot of terminal and exercises" width="800">
		</div>
	</div>
</div>

<Listings items={$tutorials} />

<Listings items={$playgrounds} title="Playgrounds" />

<Listings items={$explore} title="Explore" />

<!-- Footer -->
<Footer />

<style>
.bg-light {
	background-color: rgb(232, 240, 240) !important;
}
</style>
