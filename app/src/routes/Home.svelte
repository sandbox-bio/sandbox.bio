<script>
import { Modal } from "sveltestrap";
import { tutorials, explore } from "./stores/tutorials";
import Listings from "./components/Listings.svelte";

// State
let aboutIsOpen = false;  // Whether "About" modal is showing or not
</script>

<div class="container-fluid pt-3">
	<div class="row pb-0 pe-lg-0 bg-light align-items-center rounded-3 border shadow-lg">
		<div class="col-lg-6 p-3 p-lg-4 pt-lg-3">
			<h1 class="display-6 fw-bold lh-1">Interactive bioinformatics tutorials</h1>
			<br />
			<p class="lead">Learn how to use bioinformatics tools right from your browser.<br />Everything runs in a sandbox, so you can experiment all you want.</p>
			<div class="d-grid gap-2 d-md-flex justify-content-md-start mb-4 mt-5 mb-lg-3">
				<a href="/tutorials?id=bedtools-intro" class="btn btn-primary btn-lg px-4 me-md-2 fw-bold">Get started &rarr;</a>
			</div>
		</div>
		<div class="align-center col-lg-5 p-0 offset-lg-1 overflow-hidden shadow-lg">
			<img class="rounded-lg-3" src="images/cli.png" alt="Screenshot of terminal and exercises" width="800">
		</div>
	</div>
</div>

<Listings items={$tutorials} />

<Listings items={$explore} title="Explore" />

<!-- Footer -->
<footer class="pt-3 mt-5 mb-4 text-muted border-top">
	<div class="col-3">
		<h5>sandbox.bio</h5>
		<ul class="nav flex-column">
			<li class="nav-item mb-2"><a href="/tutorials" class="nav-link p-0">Tutorials</a></li>
			<li class="nav-item mb-2"><a href="/playground" class="nav-link p-0">Playground</a></li>
			<li class="nav-item mb-2"><a href="https://github.com/sandbox-bio/feedback/discussions" class="nav-link p-0">Feedback</a></li>
			<li class="nav-item mb-2">
				<button class="btn btn-link p-0 pb-1 text-decoration-none" on:click={() => aboutIsOpen = !aboutIsOpen}>About</button>
			</li>
		</ul>
	</div>
</footer>

<!-- "About" modal -->
<Modal body scrollable size="lg" header="About sandbox.bio" toggle={() => aboutIsOpen = !aboutIsOpen} isOpen={aboutIsOpen}>
	<p>Interactive tutorials for exploring bioinformatics command-line tools in a secure sandbox.</p>

	<p class="lead fw-bold mt-4 mb-1">How does it work?</p>
	<p>We compiled commonly-used bioinformatics tools to WebAssembly so that they can run in your browser (see our <a href="https://github.com/biowasm/biowasm" target="_blank">biowasm</a> project).</p>

	<p>To simulate a terminal environment, we implemented features such as piping (<code>|</code>), file redirection (<code>></code>, <code>>></code>), process substitution (<code>&lt;()</code>), conditional commands (<code>&&</code>, <code>||</code>), asynchronous commands (<code>&</code>), variables (<code>abc=123</code>), and autocomplete. The terminal UI is based on <a href="https://github.com/xtermjs/xterm.js/" target="_blank">xterm.js</a>.</p>

	<p>GNU Coreutils such as <code>ls</code>, <code>cat</code>, <code>grep</code>, <code>head</code>, <code>wc</code>, and <code>echo</code> were also implemented, though keep in mind we don't currently support most flags (to address that in the future, we'll compile coreutils to WebAssembly!).
		All the files you read and write to are temporarily stored in memory using Emscripten's <a href="https://emscripten.org/docs/api_reference/Filesystem-API.html" target="_blank">virtual file system</a>.</p>

	<p class="lead fw-bold mt-4 mb-1">How to contribute</p>
	<p>If you have feedback, ideas for new tutorials, or if one of your own bioinformatics tutorials could benefit from being interactive, please <a href="https://github.com/sandbox-bio/feedback/discussions" target="_blank">reach out</a>! Note that currently, only the C/C++ bioinformatics tools listed <a href="https://github.com/biowasm/biowasm#supported-tools" target="_blank">here</a> are supported.</p>

	<p class="lead fw-bold mt-4 mb-1">Author</p>
	<p>Built by <a href="https://www.robertaboukhalil.com/" target="_blank">Robert Aboukhalil</a>.</p>
</Modal>

<style>
.bg-light {
	background-color: rgb(232, 240, 240) !important;
}
</style>
