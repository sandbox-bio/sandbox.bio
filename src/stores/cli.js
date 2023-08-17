import { get, writable } from "svelte/store";
import { BUS_SERIAL_OUTPUT, DIR_TUTORIAL } from "$src/config";

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
	exec: (cmd, background = false, session = false) => {
		const command = `${cmd}\n`;
		const emulator = get(cli).emulator;

		if (background && !session) {
			emulator.keyboard_send_text(command);
		} else {
			// Execute a command without displaying it in the terminal. We could use `emulator.keyboard_send_text(command);`
			// and that would work for commands that affect the file system, but not commands that affect the current bash
			// session, e.g. `stty`, env variables, etc.
			if (session) {
				emulator.bus.listeners[BUS_SERIAL_OUTPUT] = [];
			}

			// Send command
			const chars = strToChars(command);
			chars.forEach((c) => emulator.bus.send("serial0-input", c));
		}

		// Bring back listeners after a short delay
		if (background && session) {
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
	mountFile: async (file) => {
		if (!(file instanceof File)) {
			const url = file;
			const blob = await fetch(url).then((d) => d.blob());
			file = blob;
			file.name = url.split("/").pop();
		}

		const buffer = await file.arrayBuffer();
		const view = new Uint8Array(buffer);
		const path = `${DIR_TUTORIAL}/${file.name}`;
		await get(cli).emulator.create_file(path, view);

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
