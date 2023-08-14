<script>
import { onMount } from "svelte";
import { Table, Modal, DropdownMenu, Dropdown, DropdownToggle, DropdownItem, Icon } from "sveltestrap";
import AnsiUp from "ansi_up";
import { watchResize } from "svelte-watch-resize";
import { FitAddon } from "xterm-addon-fit";
import { WebLinksAddon } from "xterm-addon-web-links";
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

let divXtermTerminal; // Xterm.js terminal
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
		initial_state: { url: "/v86/debian-state-base.bin.zst" },
		filesystem: { baseurl: "/v86/debian-9p-rootfs-flat/" },
		autostart: true,
		screen_dummy: true, // since we're using xterm.js, no need for "screen_container" div
		serial_container_xtermjs: divXtermTerminal,
		disable_mouse: true, // make sure we're still able to select text on the screen
		disable_speaker: true
	});

	$cli.emulator.bus.register("emulator-started", async () => {
		console.log("Console ready.");

		// Initialize variables and addons
		$cli.xterm = $cli.emulator.serial_adapter.term;
		$cli.addons = {
			serialize: new SerializeAddon(), // Used to export terminal to HTML
			fit: new FitAddon(), // Fit the terminal onto the screen
			links: new WebLinksAddon() // Turns text links into hyperlinks
		};

		// Initialize addons
		for (const addonName in $cli.addons) {
			$cli.xterm.loadAddon($cli.addons[addonName]);
		}

		// Make sure terminal takes up the entire div height-wise
		handleResize();

		// Initialize command line
		$cli.xterm.write("root@localhost:~# ");

		// Mount tutorial files
		for (const file of files) {
			await $cli.mount(`/${file}`);
		}
	});
});

// When window resizes, update terminal size
function handleResize() {
	if ($cli.addons.fit) {
		$cli.addons.fit.fit();

		// FIXME: for now, keep at 80 cols, but increase rows to match page height. Anything that isn't = 80 will give
		// behavior where a long command will not wrap to the next line but will overwrite the current line. It could
		// be because v86 needs to be updated to also change number of columns, but that hasn't worked so far.
		// 
		// Tried the following and it didn't work:
		//  - $cli.emulator.screen_adapter.set_size_text(dims.cols, dims.rows);
		//  - $cli.emulator.v86.cpu.devices.vga.set_size_text(dims.cols, dims.rows)
		//  - Comment out `SerialAdapterXtermJS.prototype.show` and do .open() after load addons
		//  - Changing `this.max_cols` in v86 `vga.js`
		const dims = $cli.addons.fit.proposeDimensions();
		$cli.xterm.resize(80, dims.rows);
	}
}

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

	// Mount files and show them on screen
	const paths = [];
	for (const file of files) {
		paths.push(await $cli.mount(file));
	}
	const pathsTxt = paths.join("\n\r# ");
	$cli.xterm.write(`\n\n\r\u001b[0;32m# Files mounted:\n\r# ${pathsTxt}\u001b[0m\n\n\r`);
	$cli.exec("");

	// Reset file selection (e.g. if select same file name again, should remount it because contents might be different)
	event.target.value = "";
}
</script>

<!-- Terminal -->
<div id="terminal" bind:this={divXtermTerminal} use:watchResize={handleResize}>
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
/* Xterm */
#terminal {
	opacity: 1;
	height: 85vh;
	max-height: 85vh;
	overflow: hidden;
}

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
