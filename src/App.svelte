<script>
import Home from "./Home.svelte";
import Tutorials from "./Tutorials.svelte";
import "bootstrap/dist/js/bootstrap.bundle";
import "bootstrap/dist/css/bootstrap.min.css";
</script>

<header class="d-flex flex-wrap justify-content-center py-3 mb-4 border-bottom">
	<a href="/" class="d-flex align-items-center ms-3 mb-3 mb-md-0 me-md-auto text-dark text-decoration-none">
		<span class="fs-4">🧬 sandbox.bio</span>
	</a>

	<ul class="nav nav-pills">
		<li class="nav-item"><a href="/tutorials" class="nav-link active">Tutorials</a></li>
		<li class="nav-item"><a href="/playground" class="nav-link">Playground</a></li>
		<li class="nav-item dropdown">
			<!-- svelte-ignore a11y-no-onchange -->
			<a class="nav-link dropdown-toggle" href="#" id="navbarDropdown" role="button" data-bs-toggle="dropdown" aria-expanded="false">
				Explore
			</a>
			<ul class="dropdown-menu" aria-labelledby="navbarDropdown">
				<li><a class="dropdown-item" href="https://alignment.sandbox.bio">Sequence alignment</a></li>
				<li><a class="dropdown-item" href="https://tsne.sandbox.bio">tSNE algorithm</a></li>
				<li><a class="dropdown-item" href="https://fastq.sandbox.bio">FASTQ QC metrics</a></li>
			</ul>
		</li>
	</ul>
</header>

<Tutorials />
