<script>
import Alert from "$src/components/Alert.svelte";
</script>

<svelte:head>
	<title>About - sandbox.bio</title>
</svelte:head>

<div class="col-8">
	<p class="lead">sandbox.bio features interactive tutorials for bioinformatics tools in a secure sandbox.</p>

	<Alert color="warning">
		This is an early preview of sandbox.bio v2, whose goal is to support more bioinformatics tools and tutorials than currently possible.
		<strong>This is currently in active development, so certain functionality may change or break temporarily.</strong>
	</Alert>

	<p class="lead fw-bold mt-4 mb-1">Tutorials</p>
	<p>
		Command-line tutorials are powered by <strong>a Linux OS that runs in your browser</strong>! (specifically, Debian 11). The sandbox comes
		preloaded with commonly used bioinformatics tools, so you can hit the ground running, with no setup required.
	</p>

	<p>
		This sandbox is powered by the <a href="https://github.com/copy/v86" target="_blank">v86 project</a>. Its main limitation is that tools that
		cannot be compiled to a 32-bit i686 architecture are not supported. See the
		<a href="https://github.com/sandbox-bio/v86/blob/master/tools/docker/debian/Dockerfile" target="_blank">sandbox.bio Dockerfile</a> for the latest list
		of supported tools.
	</p>

	<p class="lead fw-bold mt-4 mb-1">Playgrounds</p>

	<p>
		Playgrounds on sandbox.bio are powered by a lighter approach, using the <a href="https://biowasm.com" target="_blank">biowasm project</a>.
		Individual command-line tools were compiled to WebAssembly so they can run in the browser. The playgrounds use sample data from
		<a href="https://github.com/TheUpshot/chipotle">The UpShot</a>.
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
</div>
