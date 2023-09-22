import { get, writable } from "svelte/store";
import { BUS_SERIAL_APP_FILE, BUS_SERIAL_APP_OUTPUT, BUS_SERIAL_INPUT, BUS_SERIAL_OUTPUT, DIR_TUTORIAL } from "$src/config";

export const EXEC_MODE_TERMINAL = "terminal";
export const EXEC_MODE_TERMINAL_HIDDEN = "terminal-hidden";
export const EXEC_MODE_BUS = "bus";
const SANDBOX_END_MARKER = "__sandbox__";

function strToChars(str) {
	const chars = str.split("");
	return chars.map((d) => d.charCodeAt(0));
}

// Current tutorial
export const cli = writable({
	// -------------------------------------------------------------------------
	// State
	// -------------------------------------------------------------------------
	emulator: null,
	listeners: [],
	xterm: null,
	addons: {
		serialize: null,
		fit: null
	},

	// -------------------------------------------------------------------------
	// Utilities
	// -------------------------------------------------------------------------

	// Run a command on the command line
	exec: (cmd, { mode, callback } = { mode: EXEC_MODE_TERMINAL, callback: null }) => {
		let command = `${cmd}\n`;
		command = command.replace(/ \\ /g, ""); // for commands in tutorials split over multiple lines

		const chars = strToChars(command);
		const emulator = get(cli).emulator;

		// Before running the command, add an event listener if a callback is needed.
		// To support a callback we'll need to add an end marker and listen to the
		// UART1 port for that marker before we call the callback function.
		if (callback) {
			command = `(${cmd} && echo ${SANDBOX_END_MARKER}) > ${BUS_SERIAL_APP_FILE}\n`;
			let output = "";
			const listener = (byte) => {
				const char = String.fromCharCode(byte);
				output += char;
				const indexDoneMarker = output.indexOf(SANDBOX_END_MARKER);
				if (indexDoneMarker > -1) {
					callback(output.slice(0, indexDoneMarker));
					emulator.remove_listener(BUS_SERIAL_APP_OUTPUT, listener);
				}
			};
			emulator.add_listener(BUS_SERIAL_APP_OUTPUT, listener);
		}

		// Command executed in user's terminal
		if (mode === EXEC_MODE_TERMINAL) {
			chars.forEach((c) => emulator.bus.send(BUS_SERIAL_INPUT, c));
		}

		// Command executed on the background, skipping xterm, and not visible to the user
		if (mode === EXEC_MODE_BUS) {
			emulator.keyboard_send_text(command);
		}

		// Command executed in xterm, but not visible to the user. This is useful when you want to
		// run a command that affects the current bash session, but in the background: e.g. running
		// `stty` and defining env variables.
		if (mode === EXEC_MODE_TERMINAL_HIDDEN) {
			// Temporarily remove listeners so the xterm UI doesn't show the command
			emulator.bus.listeners[BUS_SERIAL_OUTPUT] = [];

			// Send command
			chars.forEach((c) => emulator.bus.send(BUS_SERIAL_INPUT, c));

			// Bring back listeners after a short delay
			setTimeout(() => {
				emulator.bus.listeners[BUS_SERIAL_OUTPUT] = get(cli).listeners;
			}, 500);
		}
	},

	// List files in a folder recursively (if path is to a file, returns [] if exists, undefined otherwise)
	ls: (path) => {
		const emulator = get(cli).emulator;

		// Loop through files in the current folder
		let result = [];
		const files = emulator.fs9p.read_dir(path);
		for (const file of files) {
			// If it's a folder, run ls on it recursively
			const filePath = `${path}/${file}`;
			const iNode = emulator.fs9p.SearchPath(filePath);
			if (emulator.fs9p.IsDirectory(iNode.id)) {
				result = result.concat(get(cli).ls(filePath));
			} else {
				result.push(filePath);
			}
		}

		return result;
	},

	// Mount a File object or URL to the file system
	mountFile: async (path, file) => {
		if (!(file instanceof File)) {
			const url = file;
			const blob = await fetch(url).then((d) => d.blob());
			file = blob;
			file.name = url.split("/").pop();
		}

		const buffer = await file.arrayBuffer();
		const view = new Uint8Array(buffer);
		await get(cli).createFile(`${DIR_TUTORIAL}/${path}`, view);

		return path;
	},

	// Read file from path as Uint8Array
	readFile: async (path) => {
		// Does file exist?
		const iNode = get(cli).emulator.fs9p.SearchPath(path);
		if (iNode.id === -1) {
			throw "File not found";
		}

		// If we know the file exists but it's empty, v86 incorrectly raises an exception
		try {
			return await get(cli).emulator.read_file(path);
		} catch (e) {
			if (e.message === "File not found") {
				return new Uint8Array(0);
			}
			throw message;
		}
	},

	// Create a file, given a path and contents (string or Uint8Array)
	createFile: async (path, contents) => {
		const emulator = get(cli).emulator;

		let buffer = contents;
		if (!(contents instanceof Uint8Array)) {
			buffer = new Uint8Array(contents.length);
			buffer.set(strToChars(contents));
		}

		// Create all sub-folders needed to store this file. Not using `cli.exec(mkdir, true)`
		// since we currently can't tell when a command is done running.
		let currPath = "";
		const pathElements = path.split("/").slice(0, -1);
		for (const element of pathElements) {
			currPath += `${element}/`;

			// If folder doesn't exist, create it
			const currINode = emulator.fs9p.SearchPath(currPath);
			if (currINode.id === -1) {
				emulator.fs9p.CreateDirectory(element, currINode.parentid);
			}
		}

		await emulator.create_file(path, buffer);
	}
});