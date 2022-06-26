<script>
import IDE from "../components/IDE.svelte";

let tool = "jq";
let command = JSON.stringify({"abc": "sdf", "def": "fsd", "ghi": { "a": 123, "b": 456 }}, null, 2);
let input;
let output;

// Language for the command IDE vs. input/output IDEs
$: langCmd = tool === "jq" ? "json" : "cpp";
$: langIO = tool === "jq" ? "json" : null;
// $: console.log(command);
</script>

<select bind:value={tool}>
	<option value="jq">jq</option>
	<option value="awk">awk</option>
	<option value="grep">grep</option>
	<option value="sed">sed</option>
</select>

<div class="row ide">
	<IDE lang={langCmd} code={command} on:update={d => command = d.detail} />
</div>

<div class="row">
	<div class="col-md-6 ide">
		<!-- <div bind:this={elIDEInput} /> -->
	</div>
	<div class="col-md-6 ide">
		<!-- <div bind:this={elIDEOutput} /> -->
	</div>
</div>

<!-- <button on:click={extract}>extract</button> -->

<!-- <button on:click={toggle}>toggle</button> -->

<style>
.ide {
	/* border: 1px solid #ccc; */
}
</style>
