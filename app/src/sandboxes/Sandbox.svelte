<script>
import IDE from "../components/IDE.svelte";

let tool = "jq";

let command = `.[] | select(.name == "sup")`;
let input = JSON.stringify({"abc": "sdf", "def": "fsd", "ghi": { "a": 123, "b": 456 }}, null, 2);
let output;

// Language for the command IDE vs. input/output IDEs
$: langCmd = tool === "jq" ? "json" : "cpp";
$: langIO = tool === "jq" ? "json" : null;

// For now, just a basic transformation
$: output = input + "====";
</script>

<select bind:value={tool}>
	<option value="jq">jq</option>
	<option value="awk">awk</option>
	<option value="grep">grep</option>
	<option value="sed">sed</option>
</select>

<div class="row ide ide-command mb-4">
	<IDE lang={langCmd} code={command} on:update={d => command = d.detail} />
</div>

<div class="row">
	<div class="col-md-6 ide">
		<IDE lang={langIO} code={input} on:update={d => input = d.detail} />
	</div>
	<div class="col-md-6 ide">
		<IDE lang={langIO} code={output} on:update={d => output = d.detail} editable={false} />
	</div>
</div>

<style>
.ide {
	font-size: 15px;  /* default = 16px */
	/* border: 1px solid #ccc; */
	/* overflow: scroll; */
	/* max-height: 85vh; */
}

.ide-command {
	max-height: 120px;
	overflow: scroll;
}
</style>
