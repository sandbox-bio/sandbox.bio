<script>
import { Button, Card } from "sveltestrap";
import { cli } from "$stores/cli";

export let command;
export let inline = false;

$: commandToRun = command
	.replace(
		/\\n/g,
		`
`
	)
	.replace(/ \\ /g, " ")
	.replaceAll("&lt;", "<")
	.replaceAll("&lcub;", "{")
	.replaceAll("&rcub;", "}")
$: commandPretty = command.replace(/ \\ /g, " \\ <br />&nbsp;&nbsp;&nbsp;").replace(/\\n/g, "<br>");
</script>

{#if inline}
	<Button size="sm" class="font-monospace bg-dark" on:click={() => $cli.exec(commandToRun)}>
		{@html commandPretty}
	</Button>
{:else}
	<div class="mb-3 font-monospace">
		<Card on:click={() => $cli.exec(commandToRun)} body inverse color="dark">
			{@html commandPretty}
		</Card>
	</div>
{/if}

<style>
div {
	cursor: pointer;
	font-size: 0.875em; /* Same size as <kbd> */
}
</style>
