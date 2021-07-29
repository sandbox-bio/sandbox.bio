<script>
import { xtermAddons } from "../terminal/stores";

export let command;

// Run command in the CLI!
function exec()
{
	if(!command)
			return;
	$xtermAddons.echo.setInput(command);
	$xtermAddons.echo.handleData("\r");
};
</script>

<kbd on:click={exec}>
	{command}
</kbd>

<style>
kbd:hover {
	cursor: pointer;
}
</style>
