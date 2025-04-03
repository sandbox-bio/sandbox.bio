<script>
import { goto } from "$app/navigation";
import { page } from "$app/stores";
import { onMount } from "svelte";
import { Spinner } from "sveltestrap";

// This fixes the issue where you get redirected to a blank page after logging in
onMount(() => {
	const url = $page.url.searchParams.get("url");
	const redirectTo = !url || url.startsWith("/redirect") ? "/" : url;
	goto(redirectTo);
});
</script>

<div class="text-center">
	<Spinner color="primary" size="lg" />
</div>
