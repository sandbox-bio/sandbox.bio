<script>
import { Modal } from "sveltestrap";
import Link from "$components/Link.svelte";

let aboutIsOpen = false; // Whether "About" modal is showing or not
</script>

<footer class="pt-3 mt-5 mb-4 text-muted border-top">
	<div class="col-3">
		<h5>sandbox.bio</h5>
		<ul class="nav flex-column">
			<li class="nav-item mb-2"><a href="https://github.com/sandbox-bio/sandbox.bio/discussions" class="nav-link p-0">Feedback</a></li>
			<li class="nav-item mb-2">
				<button class="btn btn-link p-0 pb-1 text-decoration-none" on:click={() => (aboutIsOpen = !aboutIsOpen)}>About</button>
			</li>
		</ul>
	</div>
</footer>

<!-- "About" modal -->
<Modal body scrollable size="lg" header="About sandbox.bio" toggle={() => (aboutIsOpen = !aboutIsOpen)} isOpen={aboutIsOpen}>
	<p>Interactive tutorials for exploring bioinformatics command-line tools in a secure sandbox.</p>

	<p class="lead fw-bold mt-4 mb-1">How does it work?</p>
	<p>
		We compiled commonly-used bioinformatics tools (and general tools such as <code>awk</code>, <code>grep</code> and <code>coreutils</code>!) to
		WebAssembly so that they can run in your browser (see our
		<a href="https://github.com/biowasm/biowasm" target="_blank" rel="noreferrer">biowasm</a> project).
	</p>

	<p>
		To simulate a terminal environment, command-line input is parsed into an abstract syntax tree and supports features such as piping (<code>|</code
		>), file redirection (<code>></code>, <code>>></code>), process substitution (<code>&lt;()</code>), conditional commands (<code>&&</code>,
		<code>||</code>), asynchronous commands (<code>&</code>), variables (<code>abc=123</code>), and autocomplete. The terminal UI is based on
		<a href="https://github.com/xtermjs/xterm.js/" target="_blank">xterm.js</a>.
	</p>

	<p>
		All the files you read and write to are temporarily stored in memory using Emscripten's <a
			href="https://emscripten.org/docs/api_reference/Filesystem-API.html"
			target="_blank"
			rel="noreferrer">virtual file system</a
		>.
	</p>

	<p>
		The <Link href="/playgrounds/awk">playgrounds</Link> use sample data from <Link href="https://github.com/TheUpshot/chipotle">The UpShot</Link>.
	</p>

	<p class="lead fw-bold mt-4 mb-1">How to contribute</p>
	<p>
		If you have feedback, ideas for new tutorials, or if one of your own bioinformatics tutorials could benefit from being interactive, please <a
			href="https://github.com/sandbox-bio/sandbox.bio/discussions"
			target="_blank"
			rel="noreferrer">reach out</a
		>!
	</p>

	<p class="lead fw-bold mt-4 mb-1">Author</p>
	<p>Built by <a href="https://robert.bio" target="_blank" rel="noreferrer">Robert Aboukhalil</a>.</p>
</Modal>
