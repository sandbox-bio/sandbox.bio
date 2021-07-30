<script>
import { Modal } from "sveltestrap";

// Props
export let alt = "";
export let src = "";

// State
let isOpen = false;  // Whether modal is showing or not
let toggle = () => isOpen = !isOpen;
</script>

<Modal body size="xl" {toggle} {isOpen}>
	<div class="text-center" style="width:100%">
		<img alt={alt} src={src} class="rounded" style="max-height:80vh; width:auto">
	</div>
</Modal>

<p>
	<img alt={alt} src={src} on:click={toggle} style="max-height:350px; width:auto">
</p>

<style>
img:hover {
	cursor: pointer;
}
</style>
