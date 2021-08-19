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
		<img class="img-fluid rounded" alt={alt} src={src}>
	</div>
</Modal>

<p class="text-center">
	<img class="img-fluid rounded img-thumbnail" alt={alt} src={src} on:click={toggle}>
</p>

<style>
img:hover {
	cursor: pointer;
}
</style>
