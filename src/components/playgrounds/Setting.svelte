<script>
import { Tooltip } from "sveltestrap";

export let tooltip;
let divSettingInteractive;
</script>

<div bind:this={divSettingInteractive} class="pt-1">
	<slot />
</div>

{#if tooltip}
	<Tooltip target={divSettingInteractive}>
		{tooltip}
	</Tooltip>
{/if}
