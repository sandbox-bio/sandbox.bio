<svelte:head>
	<title>Tutorials - sandbox.bio</title>
</svelte:head>

<slot />

<style>
:global(.language-bash) {
	white-space: pre-wrap;
}
:global(pre:has(.language-bash)) {
	border: 1px solid #ccc;
	padding: 8px;
}
</style>
