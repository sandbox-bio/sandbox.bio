<script>
import { onMount } from "svelte";
import { page } from "$app/stores";
import { goto } from "$app/navigation";
import { categories } from "$stores/tutorials";
import TutorialList from "$src/components/TutorialList.svelte";

onMount(() => {
	// Redirect URLs from sandbox.bio v1 (/tutorials?id=tutorialName&step=123 --> /tutorials/tutorialName/123)
	const oldTutorialId = $page.url.searchParams.get("id");
	const oldTutorialStep = $page.url.searchParams.get("step");

	if (oldTutorialId && oldTutorialStep) {
		goto(`/tutorials/${oldTutorialId}/${oldTutorialStep}`);
	} else if (oldTutorialId && !oldTutorialStep) {
		goto(`/tutorials/${oldTutorialId}`);
	}
});
</script>

<svelte:head>
	<title>Tutorials - sandbox.bio</title>
</svelte:head>

<h3 class="mb-3">Tutorials</h3>

<TutorialList categories={$categories} />
