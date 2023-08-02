<script>
import { onMount } from "svelte";
import { Table, Modal } from "sveltestrap";
import AnsiUp from "ansi_up";
import { V86Starter } from "$thirdparty/v86/libv86";
import { cli } from "$stores/cli";
import "xterm/css/xterm.css";

// =============================================================================
// State
// =============================================================================

export let intro = ""; // Intro string to display on Terminal once ready (optional)
export let init = ""; // Command to run to initialize the environment (optional)
export let files = []; // Files to preload on the filesystem
export let tools = TOOLS_DEFAULT; // Aioli tools to load
export let pwd = ""; // Path relative to /shared/data where user should start at

let emulator;
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
	$cli.v86 = new V86Starter({
		wasm_path: "/v86/v86.wasm",
		memory_size: 512 * 1024 * 1024,
		vga_memory_size: 8 * 1024 * 1024,
		screen_container: divScreenContainer,
		serial_container_xtermjs: divXtermTerminal,
		initial_state: { url: "/v86/debian-state-base.bin.zst" },
		filesystem: { baseurl: "/v86/debian-9p-rootfs-flat/" },
		autostart: true
	});

	$cli.v86.bus.register("emulator-started", () => {
		console.log("Console ready.");
	});
});

// =============================================================================
// Sidebar operations
// =============================================================================

// Export ANSI to HTML and open in new tab
function exportTerminal() {
	// const terminalRaw = $xtermAddons.serialize.serialize().replace(/\[[0-9]C/g, "\t");
	// const terminalHTML = "<pre>" + new AnsiUp().ansi_to_html(terminalRaw) + "</pre>";
	// const blob = new Blob([terminalHTML], { type: "text/html" });
	// const url = URL.createObjectURL(blob);
	// window.open(url);
}

// Mount local file to virtual file system
async function mountLocalFile(event) {
// 	const files = event.target.files;
// 	if (!files) {
// 		console.warn("No file specified.");
// 		return;
// 	}

// 	// Note that files that already exist will be overwritten!
// 	const paths = await $CLI.utils.mount(files);
// 	const pathsTxt = paths
// 		.map((path, i) => {
// 			if (!files[i]) return;
// 			const extra = files[i].size <= MAX_FILE_SIZE_TO_CACHE ? "" : "   (large file; won't persist on page refresh)";
// 			return `#   ${path}${extra}`;
// 		})
// 		.join("\n");
// 	input(`\n\u001b[0;32m# Files mounted:\n${pathsTxt}\u001b[0m\n\n`);

// 	// Reset file selection (e.g. if select same file name again, should remount it because contents might be different)
// 	event.target.value = "";
}
</script>

<!-- Terminal -->
<div bind:this={divXtermTerminal} style="opacity: 1; height:85vh; max-height:85vh; overflow:hidden">
	<!-- Hamburger menu for settings -->
	<div class="cli-options text-muted">
		<button class="btn btn-outline-secondary p-0 m-0 border-0" type="button" data-bs-toggle="dropdown" aria-expanded="false" on:click={() => {
			alert("sup")
		}}>
			<i class="bi bi-three-dots-vertical" />
		</button>
		<ul class="dropdown-menu">
			<li><button class="dropdown-item" on:click={() => fileInputSingle.click()}>Mount local files</button></li>
			<li><button class="dropdown-item" on:click={() => fileInputFolder.click()}>Mount local folder</button></li>
			<li><button class="dropdown-item" on:click={exportTerminal}>Export as HTML</button></li>
			<li><button class="dropdown-item" on:click={modalKbdToggle}>Keyboard Shortcuts</button></li>
		</ul>
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
	position: absolute;
	right: 50px;
	z-index: 100;
	cursor: pointer;
}

.cli-options:hover {
	color: white !important;
}
</style>
