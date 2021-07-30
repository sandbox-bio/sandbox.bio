<script>
import { xtermAddons } from "../terminal/stores";

export let command;

// Run command in the CLI!
function exec()
{
	if(!command)
		return;
	$xtermAddons.echo.setInput(command.replace(/\\/g, ""));
	$xtermAddons.echo.handleData("\r");
};
</script>

<p>
	<kbd on:click={exec}>
		{@html command.replace(/\\/g, "\\ <br />&nbsp;&nbsp;&nbsp;")}
	</kbd>
</p>

<style>
kbd:hover {
	cursor: pointer;
}
</style>
