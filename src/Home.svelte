<script>
import { config } from "../tutorials/config";
</script>

<main role="main" class="ps-4 pe-4">
	<div class="container-fluid pt-3">
		<div class="row p-4 pb-0 pe-lg-0 pt-lg-4 bg-light align-items-center rounded-3 border shadow-lg">
			<div class="col-lg-7 p-3 p-lg-4 pt-lg-3">
				<h1 class="display-6 fw-bold lh-1">Interactive bioinformatics tutorials</h1>
				<br />
				<p class="lead">Learn how to use bioinformatics tools right from your browser.<br />Everything runs in a sandbox, so you can experiment all you want.</p>
				<div class="d-grid gap-2 d-md-flex justify-content-md-start mb-4 mt-5 mb-lg-3">
					<a href="/tutorials" class="btn btn-primary btn-lg px-4 me-md-2 fw-bold">Get started &rarr;</a>
				</div>
			</div>
			<div class="col-lg-4 offset-lg-1 p-0 overflow-hidden shadow-lg">
				<img class="rounded-lg-3" src="images/cli.png" alt="" width="520">
			</div>
		</div>
	</div>

	<div class="row mt-5">
		<h3 class="pb-2">Tutorials</h3>
	</div>
	<div class="row align-items-md-stretch">
		{#each config.tutorials as tutorial}
			<div class="col-md-3">
				<div class="h-100 p-5 border rounded-3">
					<h4>{tutorial.name}</h4>
					<p>{@html tutorial.description}</p>
					<a href="/tutorials?id={tutorial.id}" class="btn btn-outline-dark text-center">Launch</a>
				</div>
			</div>
		{/each}
	</div>

	<footer class="pt-3 mt-4 text-muted border-top">
		Built by <a href="https://robertaboukhalil.com">Robert Aboukhalil</a>. Reach out if you have tutorial ideas.
	</footer>
</main>

<style>
.bg-light {
	background-color: rgb(232, 240, 240) !important;
}
</style>
