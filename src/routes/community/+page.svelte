<script>
import { Card, CardBody, CardHeader } from "sveltestrap";
</script>

<svelte:head>
	<title>Community - sandbox.bio</title>
</svelte:head>

<h3 class="my-3">Community Tutorials</h3>

<div class="row">
	<div class="col-md-6 col-lg-4">
		<div class="card h-100">
			<div class="card-header bg-primary-subtle">
				<span class="fs-5">
					<a href="/community/linux-basics" class="fw-medium text-decoration-none stretched-link">Linux Basics</a>
					(3 tutorials)
				</span>
				<p class="card-text text-muted">By the French Institute of Bioinformatics</p>
			</div>
			<div class="card-body d-flex align-items-center justify-content-center">
				<img src="/logo-ifb.svg" alt="IFB logo" height="100" />
			</div>
		</div>
	</div>
	<div class="col-md-6 col-lg-4">
		<div class="card h-100">
			<div class="card-header bg-primary-subtle">
				<span class="fs-5">
					<a href="/community/data-processing" class="fw-medium text-decoration-none stretched-link">Data and Text Processing</a>
					(8 tutorials)
				</span>
				<p class="card-text">By Francisco M. Couto</p>
			</div>
			<div class="card-body d-flex align-items-center justify-content-center">
				<img src="/logo-book-data-processing.jpg" alt="Book cover" width="150" />
			</div>
		</div>
	</div>
</div>

<style>
.card:hover {
	border: 1px solid var(--bs-primary);
}
</style>
