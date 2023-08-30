<script>
import { onMount } from "svelte";
import { Table, Modal, DropdownMenu, Dropdown, DropdownToggle, DropdownItem, Icon, Spinner } from "sveltestrap";
import AnsiUp from "ansi_up";
import { watchResize } from "svelte-watch-resize";
import { FitAddon } from "xterm-addon-fit";
import { WebLinksAddon } from "xterm-addon-web-links";
import { SerializeAddon } from "xterm-addon-serialize";
import { V86Starter } from "$thirdparty/v86/libv86";
import { EXEC_MODE_TERMINAL, EXEC_MODE_TERMINAL_HIDDEN, cli } from "$stores/cli";
import { tutorial } from "$stores/tutorial";
import { LocalState, log } from "$src/utils";
import { BUS_SERIAL_OUTPUT, DIR_TUTORIAL, DIR_TUTORIAL_SHORT, LOGGING_DEBUG, LOGGING_INFO, MAX_FILE_SIZE_TO_CACHE } from "$src/config";
import "xterm/css/xterm.css";

// =============================================================================
// State
// =============================================================================

export let files = []; // Files to preload on the filesystem
export let intro = ""; // Intro string to display on Terminal once ready (optional)
export let init = ""; // Command to run to initialize the environment (optional)
export let tools; // Aioli tools to load

let loading = false;
let divXtermTerminal; // Xterm.js terminal
let inputMountFiles; // Hidden HTML file input element for mounting local file
let inputMountFolder; // Hidden HTML file input element for mounting local folder
let modalKbdOpen = false; // Set to true when the shortcuts modal is open
let modalKbdToggle = () => (modalKbdOpen = !modalKbdOpen);
let timerSyncFS; // JS timeout used to sync filesystem contents

// =============================================================================
// Initialization
// =============================================================================

onMount(initialize);

function initialize() {
	console.log("Initializing terminal...");
	loading = true;

	if (timerSyncFS) {
		clearTimeout(timerSyncFS);
	}

	$cli.emulator = new V86Starter({
		wasm_path: "https://assets.sandbox.bio/v86/v86.wasm",
		memory_size: 512 * 1024 * 1024,
		initial_state: { url: "https://assets.sandbox.bio/v86/debian-state-base.bin.zst" },
		filesystem: { baseurl: "https://assets.sandbox.bio/v86/debian-9p-rootfs-flat/" },
		autostart: true,
		screen_dummy: true, // since we're using xterm.js, no need for "screen_container" div
		serial_container_xtermjs: divXtermTerminal,
		disable_mouse: true, // make sure we're still able to select text on the screen
		disable_speaker: true,
		uart1: true // we'll use serial port 1 to communicate with v86 from JavaScript
	});

	$cli.emulator.bus.register("emulator-loaded", async () => {
		$cli.xterm = $cli.emulator.serial_adapter.term;
		$cli.listeners = $cli.emulator.bus.listeners[BUS_SERIAL_OUTPUT];

		// Make sure everything loaded correctly. If not, try again.
		// Otherwise, get issues where `term` variable is null and waiting for it to be set does not help.
		if (!$cli.xterm) {
			initialize();
			return;
		}
		console.log("Terminal ready.", $cli);

		// Initialize addons
		$cli.addons = {
			serialize: new SerializeAddon(), // Used to export terminal to HTML
			fit: new FitAddon(), // Fit the terminal onto the screen
			links: new WebLinksAddon() // Turns text links into hyperlinks
		};
		for (const addonName in $cli.addons) {
			$cli.xterm.loadAddon($cli.addons[addonName]);
		}
		// Make sure terminal takes up the entire div height-wise
		handleResize(true);

		// Initialize command line
		$cli.xterm.write(`root@localhost:${DIR_TUTORIAL_SHORT}# `);

		// Mount tutorial files
		await mountTutorialFiles();
		// Mount previously synced FS (user's FS takes precedence over tutorial files)
		await fsLoad();
		// Start syncing FS
		fsSync();

		loading = false;
	});
}

// When window resizes, update terminal size
function handleResize(hidden = false) {
	if ($cli.addons.fit) {
		$cli.addons.fit.fit();

		// If we resize the terminal's number of rows/cols on xterm.js, we also need to update those
		// values for the actual terminal itself. Otherwise, the following issues arise:
		// - Long commands don't wrap to the next line and start overwriting the start of the command
		// - Editing previously run long-commands shows odd spacing behavior
		// - TUIs like `top` and `vim` don't load in full screen
		const dims = $cli.addons.fit.proposeDimensions();
		log(LOGGING_INFO, "Resize terminal", dims);
		$cli.exec(`stty rows ${dims.rows} cols ${dims.cols}`, {
			mode: hidden ? EXEC_MODE_TERMINAL_HIDDEN : EXEC_MODE_TERMINAL
		});
	}
}

// =============================================================================
// File system sync
// =============================================================================

async function fsSync() {
	// If navigate away from tutorial, stop syncing
	if (!$tutorial.id) return;
	log(LOGGING_DEBUG, "Saving FS state...");

	await fsSave();
	timerSyncFS = setTimeout(fsSync, 1000);
}

// Save FS state to localforage
async function fsSave() {
	const fs = [];
	const paths = $cli.ls(DIR_TUTORIAL);
	for (const path of paths) {
		const contents = await $cli.readFile(path);
		if (contents.length <= MAX_FILE_SIZE_TO_CACHE) {
			fs.push({ path, contents });
		}
	}
	await LocalState.setFS($tutorial.id, fs);
}

// Load FS state from localforage
async function fsLoad() {
	const files = await LocalState.getFS($tutorial.id);
	for (const file of files) {
		await $cli.createFile(file.path, file.contents);
	}
}

// =============================================================================
// Sidebar operations
// =============================================================================

// Mount tutorial files
async function mountTutorialFiles() {
	for (const file of files) {
		await $cli.mountFile(file, `/data/${$tutorial.id}/${file}`)
	}
}

// Export ANSI to HTML and open in new tab
function exportHTML() {
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
		paths.push(await $cli.mountFile(file.name, file));
	}
	const pathsTxt = paths.join("\n\r# ");
	$cli.xterm.write(`\n\n\r\u001b[0;32m# Files mounted:\n\r# ${pathsTxt}\u001b[0m\n\n\r`);
	$cli.exec("");

	// Reset file selection (e.g. if select same file name again, should remount it because contents might be different)
	event.target.value = "";
}
</script>

<!-- Terminal -->
<div id="terminal" bind:this={divXtermTerminal} use:watchResize={handleResize} class:opacity-25={loading}>
	{#if loading}
		<Spinner color="light" type="border" style="position:absolute" />
	{/if}
	<div class="cli-options">
		<Dropdown autoClose={true}>
			<DropdownToggle color="dark" size="sm">
				<Icon name="three-dots-vertical" />
			</DropdownToggle>
			<DropdownMenu>
				<DropdownItem on:click={mountTutorialFiles}>Reset tutorial files</DropdownItem>
				<DropdownItem on:click={() => inputMountFiles.click()}>Mount local files</DropdownItem>
				<DropdownItem on:click={() => inputMountFolder.click()}>Mount local folder</DropdownItem>
				<DropdownItem on:click={exportHTML}>Export as HTML</DropdownItem>
				<DropdownItem on:click={modalKbdToggle}>Keyboard Shortcuts</DropdownItem>
			</DropdownMenu>
		</Dropdown>
	</div>
</div>

<!-- Hidden input file for mounting local files -->
<input type="file" on:change={mountLocalFile} bind:this={inputMountFiles} style="display:none" multiple />
<input type="file" on:change={mountLocalFile} bind:this={inputMountFolder} style="display:none" multiple webkitdirectory />

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
