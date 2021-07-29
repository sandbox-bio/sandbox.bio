<script>
import Terminal from "./terminal/Terminal.svelte";
</script>

<div class="container-fluid pb-3">
	<div class="d-grid gap-3" style="grid-template-columns: 1fr 3fr;">
		<div class="bg-light border rounded-3 p-2">
			
		</div>
		<div id="terminal-wrapper" class="border rounded-3 p-2">
			<Terminal />
		</div>
	</div>
</div>

<style>
#terminal-wrapper {
	background-color: black;
}
</style>
