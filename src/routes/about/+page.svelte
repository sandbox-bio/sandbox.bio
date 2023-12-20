<svelte:head>
	<title>About - sandbox.bio</title>
</svelte:head>

<div class="col-8">
	<p class="lead">sandbox.bio features interactive tutorials for bioinformatics tools in a secure sandbox.</p>

	<p class="lead fw-bold mt-4 mb-1">Tutorials</p>
	<p>
		Command-line tutorials feature a <strong>Linux operating system running in your browser!</strong> The sandbox runs Debian 12 (Bookworm), and comes
		preloaded with commonly used bioinformatics tools, so you can hit the ground running, with no setup required. To use sandbox.bio, you need a computer
		with >1GB of RAM and a connection to the internet.
	</p>

	<p>
		The Debian environment is powered by the <a href="https://github.com/copy/v86" target="_blank">v86 project</a>. For a tool to be available in this
		environment, it must be compiled to a 32-bit i686 architecture, and only use supported SIMD instructions (SSE, SSE2, SSE3). Currently, JVM-based
		tools are not supported because they require large downloads, and are not performant in this environment. See the
		<a href="https://github.com/sandbox-bio/v86/blob/master/tools/docker/debian/Dockerfile" target="_blank">sandbox.bio Dockerfile</a> for the latest list
		of available tools. Please reach out if you're interested in adding more tools to the sandbox.
	</p>

	<p>
		The file system available in the sandbox is a virtual file system, and is not connected to your computer's files. Each tutorial gets its own <code
			>~/tutorial</code
		> folder, whose contents is saved to your browser's storage every few seconds, so refreshing the page restores the previous state of the file system.
		This browser storage is temporary, so do not store important files in this environment.
	</p>

	<p class="lead fw-bold mt-4 mb-1">Playgrounds</p>

	<p>
		Playgrounds on sandbox.bio are powered by <a href="https://biowasm.com" target="_blank">biowasm</a>. Individual command-line tools were compiled
		to WebAssembly so they can run in the browser. The playgrounds use sample data from
		<a href="https://github.com/TheUpshot/chipotle">The UpShot</a>.
	</p>

	<p class="lead fw-bold mt-4 mb-1">How to contribute</p>
	<p>
		If you have feedback, please reach out on <a href="https://github.com/sandbox-bio/sandbox.bio/discussions/categories/general">GitHub Discussions</a>. If you're interested in contributing tutorials, please refer to our <a href="https://github.com/sandbox-bio/sandbox.bio/blob/main/CONTRIBUTING.md">Contributing documentation</a>.
	</p>

	<p class="lead fw-bold mt-4 mb-1">Author</p>
	<p>Built by <a href="https://robert.bio" target="_blank" rel="noreferrer">Robert Aboukhalil</a>.</p>
</div>
