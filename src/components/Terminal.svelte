<script>
import { onDestroy, onMount } from "svelte";
import { Table, Modal, DropdownMenu, Dropdown, DropdownToggle, DropdownItem, Icon, Spinner } from "sveltestrap";
import AnsiUp from "ansi_up";
import { watchResize } from "svelte-watch-resize";
import { FitAddon } from "xterm-addon-fit";
import { WebLinksAddon } from "xterm-addon-web-links";
import { SerializeAddon } from "xterm-addon-serialize";
import { env } from "$env/dynamic/public";
import { V86 } from "$thirdparty/v86/libv86";
import { EXEC_MODE_BUS, EXEC_MODE_TERMINAL_HIDDEN, cli } from "$stores/cli";
import { LocalState, log } from "$src/utils";
import {
	BUS_OUTPUT_CUSTOM_COMMAND,
	BUS_INPUT,
	BUS_OUTPUT,
	DEBIAN_STATE_ID,
	DIR_TUTORIAL,
	LOGGING_INFO,
	MAX_FILE_SIZE_TO_CACHE,
	URL_ASSETS
} from "$src/config";
import "xterm/css/xterm.css";
const DEBUG = false;
const SYNC_FS = true;

// =============================================================================
// State
// =============================================================================

export let terminalId = "terminal";
export let files = []; // Files to preload on the FS from /data/<tutorial>
export let assets = []; // Files to preload on the FS from assets.sandbox.bio/tutorials/<tutorial>
export let init = ""; // Command to run to initialize the environment (optional)
export let tools = []; // For these tools, pre-download .bin files (optional)
export let intro = ""; // Intro string to display on Terminal once ready (optional; not currently used in any tutorial)

let loading = true; // Loading the terminal
let loadingStatus = []; // Loading progress to show (each element = 1 line)
let mounted = false; // Component is mounted and ready to go
let divXtermTerminal; // Xterm.js terminal
let inputMountFiles; // Hidden HTML file input element for mounting local file
let inputMountFolder; // Hidden HTML file input element for mounting local folder
let modalKbdOpen = false; // Set to true when the shortcuts modal is open
let modalKbdToggle = () => (modalKbdOpen = !modalKbdOpen);
let timerSyncFS; // JS timeout used to sync filesystem contents
let timerWaitForPrompt; // Wait for root@localhost prompt to be visible

// Environments where sandbox.bio is deployed
const environments = {
	localhost: {
		url: "",
		v86: ""
	},
	"stg.sandbox.bio": {
		url: URL_ASSETS,
		v86: "stg/"
	},
	"beta.sandbox.bio": {
		url: URL_ASSETS,
		v86: "prd/"
	},
	"sandbox.bio": {
		url: URL_ASSETS,
		v86: "prd/"
	}
};

function getEnvironmentInfo() {
	// Uncomment for debugging assets on other environments
	// return environments["stg.sandbox.bio"];

	// If running tests on GitHub or using preview branches, run them on prd assets, despite being on localhost
	if (env.PUBLIC_USE_PRD_ASSETS === "true" || window.location.hostname.endsWith(".sandbox-bio.pages.dev")) return environments["sandbox.bio"];

	// Otherwise, use hostname
	const envInfo = environments[window.location.hostname];
	if (!envInfo) throw "Unrecognized deploy environment";
	return envInfo;
}

// =============================================================================
// Initialization
// =============================================================================

// Needs to be mounted or get errors on first mount
$: if (mounted && terminalId) initialize(terminalId);
onMount(() => (mounted = true));
onDestroy(cleanupTimers);

function cleanupTimers() {
	clearTimeout(timerSyncFS);
	clearInterval(timerWaitForPrompt);
}

function addLoadingStatus(msg) {
	loadingStatus = [...loadingStatus, msg];
}

function initialize(id) {
	console.log("Initializing terminal...", id);
	loading = true;
	console.time("initialize");
	addLoadingStatus("Setting up your terminal...");

	// Cleanup
	cleanupTimers();
	if ($cli.emulator) $cli.emulator.destroy();

	// Create emulator
	const envInfo = getEnvironmentInfo();
	$cli.emulator = new V86({
		wasm_path: `/v86/v86.wasm`,
		memory_size: 1024 * 1024 * 1024,
		initial_state: { url: `${envInfo.url}/v86/${envInfo.v86}debian-state-${DEBIAN_STATE_ID}.bin.zst` },
		filesystem: { baseurl: `${envInfo.url}/v86/${envInfo.v86}debian-9p-rootfs-flat/` },
		autostart: true,
		screen_dummy: DEBUG, // since we're using xterm.js, no need for "screen_container" div
		screen_container: DEBUG ? document.getElementById("screen_container") : null,
		serial_container_xtermjs: divXtermTerminal,
		disable_mouse: true, // make sure we're still able to select text on the screen
		disable_speaker: true,
		uart1: true, // see config.js for how we're using serial ports to communicate between JS and v86
		uart2: true
	});

	// Listen on the special ttyS1 port for communication from within the emulator
	let output = "";
	$cli.emulator.add_listener(BUS_OUTPUT_CUSTOM_COMMAND, async (byte) => {
		const char = String.fromCharCode(byte);
		if (char !== "\n") {
			output += char;
		} else {
			try {
				const command = JSON.parse(output);
				const params = command.params;
				console.log("Command:", command);

				// Open file contents in a new tab
				if (command.type === "open") {
					const contents = await $cli.emulator.read_file(params.path);
					const blob = new Blob([contents], { type: params.path.endsWith(".html") ? "text/html" : "text/plain" });
					const url = URL.createObjectURL(blob);
					window.open(url);

					// Launch file download
				} else if (command.type === "download") {
					const contents = await $cli.emulator.read_file(params.path);
					const blob = new Blob([contents], { type: "application/octet-stream" });
					const url = URL.createObjectURL(blob);

					// Create link element to customize the filename, otherwise it's a UUID.
					const fileLink = document.createElement("a");
					fileLink.href = url;
					fileLink.download = params.path.split("/").pop();
					fileLink.click();

					// Make fetch call and return result
					// Storing result in file; when try to store result in /dev/ttyS2, it adds a bunch of "\n" and skips the first few bytes
				} else if (command.type === "fetch" && params.url && params.output && params.sentinel) {
					// Download file and save to FS
					const buffer = await fetch(params.url).then((d) => d.arrayBuffer());
					await $cli.emulator.create_file(params.output, new Uint8Array(buffer));
					await $cli.emulator.create_file(params.sentinel, new Uint8Array([]));
				}
			} catch (e) {
				console.log("Received:", output);
				console.error(e);
			}

			output = "";
		}
	});

	// Listen for outputs
	let initial_screen = "";
	const listenerWaitForPrompt = async (byte) => (initial_screen += String.fromCharCode(byte));
	$cli.emulator.add_listener(BUS_OUTPUT, listenerWaitForPrompt);

	// Prepare terminal environment
	$cli.emulator.bus.register("emulator-loaded", async () => {
		$cli.xterm = $cli.emulator.serial_adapter.term;
		$cli.listeners = $cli.emulator.bus.listeners[BUS_OUTPUT];

		// Make sure everything loaded correctly. If not, try again.
		// Otherwise, get issues where `term` variable is null and waiting for it to be set does not help.
		if (!$cli.xterm) {
			loading = false;
			console.warn("Could not load terminal; serial_adapter not defined.");
			initialize(terminalId);
			return;
		}

		// Initialize addons
		$cli.addons = {
			serialize: new SerializeAddon(), // Used to export terminal to HTML
			fit: new FitAddon(), // Fit the terminal onto the screen
			links: new WebLinksAddon() // Turns text links into hyperlinks
		};
		for (const addonName in $cli.addons) {
			$cli.xterm.loadAddon($cli.addons[addonName]);
		}
		console.log("Terminal ready.", $cli);

		// Mount tutorial files and previously synced FS (user's FS overrides default tutorial files)
		addLoadingStatus("Loading tutorial files...");
		await mountTutorialFiles();
		await fsLoad();

		// Preload tools so by the time the user needs them, they are cached. We're done fetching
		// data from the server so won't compete with other fetch requests.  We use "&" to download
		// files in parallel as much as possible. The alternative would be to download .bin
		// files directly but we'd have to generate a list of .bin from `debian-base-fs.json`,
		// which is prone to changes.
		$cli.exec(`sync & echo & ls & ll & pwd`, { mode: EXEC_MODE_BUS });
		$cli.exec(tools.map((t) => `timeout 2 ${t}`).join(" & "), { mode: EXEC_MODE_BUS });

		// Run initialization commands
		addLoadingStatus("Initializing environment...");
		$cli.exec(init);
		// Set initial terminal size, otherwise sometimes doesn't call that function at load time
		handleResize(true);
		// Focus cursor on command line
		$cli.xterm.focus();

		// Sync date and time (otherwise continues from date/time from last boot)
		$cli.exec(`date -s "${new Date().toString()}"`, { mode: EXEC_MODE_BUS });

		// Make sure root@localhost prompt shows up on screen
		addLoadingStatus("Putting the finishing touches...");
		timerWaitForPrompt = setInterval(() => {
			if (!initial_screen.includes("root@localhost")) {
				$cli.exec("");
				// Press Ctrl + L (key code 12) to show the prompt but without extra lines above it
				$cli.emulator.bus.send(BUS_INPUT, 12);
				if (intro) $cli.exec(intro);
			} else {
				$cli.emulator.remove_listener(BUS_OUTPUT, listenerWaitForPrompt);
				clearInterval(timerWaitForPrompt);
				loading = false;
				console.timeEnd("initialize");

				// Start syncing FS
				fsSync();
			}
		}, 200);
	});
}

// When window resizes, update terminal size
let currDims = { cols: null, rows: null };
function handleResize(firstTime = false) {
	if (loading && !firstTime) return;
	if (!$cli.addons.fit) return;

	$cli.addons.fit.fit();

	// If we resize the terminal's number of rows/cols on xterm.js, we also need to update those
	// values for the actual terminal itself. Otherwise, the following issues arise:
	// - Long commands don't wrap to the next line and start overwriting the start of the command
	// - Editing previously run long-commands shows odd spacing behavior
	// - TUIs like `top` and `vim` don't load in full screen
	const dims = $cli.addons.fit.proposeDimensions();
	if (!dims?.cols || !dims?.rows || (currDims.cols === dims.cols && currDims.rows === dims.rows)) return;
	currDims = dims;

	// Limitation: this doesn't work if you're inside vim/less/etc, or halfway through a command
	// before resizing, but that should be less likely.
	if (firstTime) log(LOGGING_INFO, "Set terminal size", dims);
	else log(LOGGING_INFO, "Resize terminal", dims);
	$cli.exec(`stty rows ${dims.rows} cols ${dims.cols}`, {
		mode: firstTime ? EXEC_MODE_TERMINAL_HIDDEN : null
	});
}

// =============================================================================
// File system sync
// =============================================================================

async function fsSync() {
	if (!SYNC_FS) return;
	try {
		await fsSave();
		timerSyncFS = setTimeout(fsSync, 2000);
	} catch (error) {
		console.warn(error);
	}
}

// Save FS state to localforage
async function fsSave() {
	const id = terminalId;

	// Clear the cache before we save the file system state (otherwise, some files get stored as empty files)
	await $cli.clearCache();

	// Export FS state
	const files = $cli.ls(DIR_TUTORIAL);

	for (const file of files) {
		// If file, get contents
		if (!file.isDir) {
			const contents = await $cli.readFile(file.path);
			if (contents.length <= MAX_FILE_SIZE_TO_CACHE) {
				file.contents = contents;
			}
		}
	}

	// Only sync FS if did not switch tutorials in the middle of syncing the FS
	if (terminalId === id) await LocalState.setFS(terminalId, files);
}

// Load FS state from localforage
async function fsLoad() {
	const files = await LocalState.getFS(terminalId);
	for (const file of files) {
		if (file.isDir) {
			await $cli.createFolder(file.path);
		} else if (file.contents) {
			await $cli.createFile(file.path, file.contents);
		} else {
			console.warn(`Skipping file ${file?.path} because too large to save in browser.`);
		}
	}
}

// =============================================================================
// Sidebar operations
// =============================================================================

// Mount tutorial files
async function mountTutorialFiles() {
	// Clear cache before we mount files, otherwise they don't change!
	await $cli.clearCache();

	// Mount files stored in this repo
	for (const fileName of files) {
		const url = `/data/${terminalId}/${fileName}`;
		await $cli.mountFile(fileName, url);
	}
	// Mount files stored in assets.sandbox.bio because of their size
	for (const fileName of assets || []) {
		const url = `https://assets.sandbox.bio/tutorials/${terminalId}/${fileName}`;
		await $cli.mountFile(fileName, url);
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
		// Try to preserve folder structure
		const fileName = file?.webkitRelativePath || file.name;
		paths.push(await $cli.mountFile(fileName, file));
	}
	const pathsTxt = paths.join("\n\r# ");
	$cli.xterm.write(`\n\n\r\u001b[0;32m# Files mounted:\n\r# ${pathsTxt}\u001b[0m\n\n\r`);
	$cli.exec("");

	// Reset file selection (e.g. if select same file name again, should remount it because contents might be different)
	event.target.value = "";
}
</script>

<!-- Debug ttyS1 -->
{#if DEBUG}
	<div id="screen_container">
		<div id="screen" />
		<canvas id="vga" />
		<div style="position: absolute; top: 0; z-index: 10">
			<textarea class="phone_keyboard" />
		</div>
	</div>
{/if}

<!-- Terminal -->
{#if loading}
	<div style="position:absolute">
		<div class="text-light font-monospace small">
			<Spinner color="light" type="border" size="sm" />
			{#each loadingStatus as status, i}
				<span class:ps-4={i > 0}>
					{status}
				</span>
				{#if i < loadingStatus.length - 1}
					<span class="text-success fw-bold">ok</span>
				{/if}
				<br />
			{/each}
		</div>
	</div>
{/if}
<div id="terminal" bind:this={divXtermTerminal} use:watchResize={handleResize} class:opacity-0={loading}>
	<div class="cli-options">
		<Dropdown autoClose={true}>
			<DropdownToggle color="dark" size="sm">
				<Icon name="three-dots-vertical" />
			</DropdownToggle>
			<DropdownMenu>
				{#if terminalId !== "playground"}
					<DropdownItem on:click={mountTutorialFiles}>Reset tutorial files</DropdownItem>
				{/if}
				<DropdownItem on:click={() => inputMountFiles.click()}>Mount local files</DropdownItem>
				<DropdownItem on:click={() => inputMountFolder.click()}>Mount local folder</DropdownItem>
				<DropdownItem on:click={exportHTML}>Export as HTML</DropdownItem>
				<DropdownItem on:click={modalKbdToggle}>Shortcuts</DropdownItem>
			</DropdownMenu>
		</Dropdown>
	</div>
</div>

<!-- Hidden input file for mounting local files -->
<input type="file" on:change={mountLocalFile} bind:this={inputMountFiles} style="display:none" multiple />
<input type="file" on:change={mountLocalFile} bind:this={inputMountFolder} style="display:none" multiple webkitdirectory />

<!-- Keyboard Shortcuts Modal -->
<Modal body header="Shortcuts" isOpen={modalKbdOpen} toggle={modalKbdToggle}>
	<h4>Keyboard Shortcuts</h4>
	<Table>
		<thead>
			<tr>
				<th width="40%">Shortcut</th>
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

	<h4 class="mt-5">Special Commands</h4>
	<Table>
		<thead>
			<tr>
				<th width="40%">Command</th>
				<th>Action</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td><code>download file.txt</code></td>
				<td>Downloads <code>file.txt</code> to your computer</td>
			</tr>
			<tr>
				<td><code>open file.html</code></td>
				<td>Opens <code>file.html</code> in a new tab</td>
			</tr>
			<tr>
				<td><code>curl -O https://.../file.txt</code></td>
				<td>
					Downloads <code>file.txt</code> to the sandbox. Note that because of browser limitations, <code>curl</code> only supports URLs that are CORS-enabled.
				</td>
			</tr>
		</tbody>
	</Table>
</Modal>

<style>
/* Xterm */
#terminal {
	height: 85vh;
	max-height: 85vh;
	overflow: hidden;
}

/* Hamburger menu */
.cli-options {
	position: absolute;
	right: 0;
	z-index: 100;
}

.cli-options:hover {
	color: white !important;
}
</style>
