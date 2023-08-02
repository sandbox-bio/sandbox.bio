<script>
import { onMount } from "svelte";
import { Table, Modal, DropdownMenu, Dropdown, DropdownToggle, DropdownItem, Icon } from "sveltestrap";
import AnsiUp from "ansi_up";
import { SerializeAddon } from "xterm-addon-serialize";
import { V86Starter } from "$thirdparty/v86/libv86";
import { cli } from "$stores/cli";
import "xterm/css/xterm.css";

// =============================================================================
// State
// =============================================================================

export let intro = ""; // Intro string to display on Terminal once ready (optional)
export let init = ""; // Command to run to initialize the environment (optional)
export let files = []; // Files to preload on the filesystem
export let tools; // Aioli tools to load
export let pwd = ""; // Path relative to /shared/data where user should start at

let divScreenContainer;
let divXtermTerminal;

let fileInputSingle; // Hidden HTML file input element for mounting local file
let fileInputFolder; // Hidden HTML file input element for mounting local folder
let modalKbdOpen = false; // Set to true when the shortcuts modal is open
let modalKbdToggle = () => (modalKbdOpen = !modalKbdOpen);

// =============================================================================
// Initialization
// =============================================================================

// On mount
onMount(() => {
	$cli.emulator = new V86Starter({
		wasm_path: "/v86/v86.wasm",
		memory_size: 512 * 1024 * 1024,
		vga_memory_size: 8 * 1024 * 1024,
		screen_container: divScreenContainer,
		serial_container_xtermjs: divXtermTerminal,
		initial_state: { url: "/v86/debian-state-base.bin.zst" },
		filesystem: { baseurl: "/v86/debian-9p-rootfs-flat/" },
		autostart: true
	});

	$cli.emulator.bus.register("emulator-started", () => {
		console.log("Console ready.");

		// Initialize variables and addons
		$cli.xterm = $cli.emulator.serial_adapter.term;
		$cli.addons.serialize = new SerializeAddon();
		$cli.xterm.loadAddon($cli.addons.serialize);
	});
});

// =============================================================================
// Sidebar operations
// =============================================================================

// Export ANSI to HTML and open in new tab
function exportTerminal() {
	const terminalRaw = $cli.addons.serialize.serialize();
	const terminalHTML = "<pre>" + new AnsiUp().ansi_to_html(terminalRaw) + "</pre>";
	const blob = new Blob([terminalHTML], { type: "text/html" });
	const url = URL.createObjectURL(blob);
	window.open(url);
}

// Mount local file to virtual file system
async function mountLocalFile(event) {
	const files = event.target.files;
	if (!files) {
		console.warn("No file specified.");
		return;
	}

	// Mount files (synchronously; couldn't get AsyncFileBuffer working)
	for(const file of files) {
		var loader = new $cli.emulator.v86util.SyncFileBuffer(file);
		loader.onload = function() {
			loader.get_buffer(async function(buffer) {
				await $cli.emulator.create_file(`/root/${file.name}`, new Uint8Array(buffer));
			});
		};
		loader.load();
	}

	// Reset file selection (e.g. if select same file name again, should remount it because contents might be different)
	event.target.value = "";
}
</script>

<!-- Terminal -->
<div bind:this={divXtermTerminal}>
	<div class="cli-options">
		<Dropdown autoClose={true}>
			<DropdownToggle color="dark" size="sm">
				<Icon name="three-dots-vertical" />
			</DropdownToggle>
			<DropdownMenu>
				<DropdownItem on:click={() => fileInputSingle.click()}>Mount local files</DropdownItem>
				<DropdownItem on:click={() => fileInputFolder.click()}>Mount local folder</DropdownItem>
				<DropdownItem on:click={exportTerminal}>Export as HTML</DropdownItem>
				<DropdownItem on:click={modalKbdToggle}>Keyboard Shortcuts</DropdownItem>
			</DropdownMenu>
		</Dropdown>
	</div>
</div>

<!-- Hidden v86 terminal -->
<div style="display: none" bind:this={divScreenContainer}>
	<div id="screen" />
	<canvas id="vga" />
	<div style="position: absolute; top: 0; z-index: 10">
		<textarea class="phone_keyboard" />
	</div>
</div>

<!-- Hidden input file for mounting local files -->
<input type="file" on:change={mountLocalFile} bind:this={fileInputSingle} style="display:none" multiple />
<input type="file" on:change={mountLocalFile} bind:this={fileInputFolder} style="display:none" multiple webkitdirectory />

<!-- Keyboard Shortcuts Modal -->
<Modal body header="Keyboard Shortcuts" isOpen={modalKbdOpen} toggle={modalKbdToggle}>
	<Table>
		<thead>
			<tr>
				<th>Shortcut</th>
				<th>Action</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td><code>Ctrl + L</code></td>
				<td>Clear</td>
			</tr>
			<tr>
				<td><code>Ctrl + A</code></td>
				<td>Go to start of line</td>
			</tr>
			<tr>
				<td><code>Ctrl + E</code></td>
				<td>Go to end of line</td>
			</tr>
			<tr>
				<td><code>Ctrl + W</code></td>
				<td>Delete previous word</td>
			</tr>
			<tr>
				<td><code>Alt + Left</code></td>
				<td>Go to previous word</td>
			</tr>
			<tr>
				<td><code>Alt + Right</code></td>
				<td>Go to following word</td>
			</tr>
		</tbody>
	</Table>
</Modal>

<style>
/* Hamburger menu */
.cli-options {
	position: relative;
	float: right;
	z-index: 100;
}

.cli-options:hover {
	color: white !important;
}
</style>
