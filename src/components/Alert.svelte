<script>
import { Alert } from "sveltestrap";

export let color = "info";
</script>

<Alert {color} fade={false}>
	<slot />
</Alert>
