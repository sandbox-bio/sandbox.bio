<script>
import { Card } from "sveltestrap";
import { xtermAddons } from "terminal/xterm";

export let command;
export let inline = false;

// Run command in the CLI!
function exec() {
	if(!command)
		return;
	$xtermAddons.echo.handleData(command.replace(/ \\ /g, " "))
	$xtermAddons.echo.handleData("\r");
};
</script>

{#if inline}
	<kbd on:click={exec}>
		{@html command.replace(/ \\ /g, " \\ <br />&nbsp;&nbsp;&nbsp;")}
	</kbd>

{:else}
	<div class="cursor-pointer mb-3 font-monospace">
		<Card on:click={exec} body inverse color="dark">
			{@html command.replace(/ \\ /g, " \\ <br />&nbsp;&nbsp;&nbsp;")}
		</Card>
	</div>
{/if}

<style>
kbd, div {
	cursor: pointer;
}

div {
	font-size: .875em;  /* Same size as <kbd> */
}
</style>
