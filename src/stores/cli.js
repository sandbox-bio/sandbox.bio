import { get, writable } from "svelte/store";

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
	xterm: null,
	addons: {
		serialize: null,
		fit: null
	},

	// -------------------------------------------------------------------------
	// Utilities
	// -------------------------------------------------------------------------

	// Run a command on the command line
	exec: (cmd) => {
		const v86 = get(cli).emulator;
		const chars = strToChars(`${cmd}\n`);
		chars.forEach((c) => v86.bus.send("serial0-input", c));
	},

	// Mount a File object or URL to the file system
	mount: async (file, folder = "/root") => {
		if (!(file instanceof File)) {
			const url = file;
			const blob = await fetch(url).then((d) => d.blob());
			file = blob;
			file.name = url.split("/").pop();
		}

		const buffer = await file.arrayBuffer();
		const view = new Uint8Array(buffer);
		const path = `${folder}/${file.name}`;
		await get(cli).emulator.create_file(path, view);

		return path;
	},

	// -------------------------------------------------------------------------
	// Not yet used
	// -------------------------------------------------------------------------

	// List files in a folder
	ls: (path) => {
		return get(cli).emulator.fs9p.read_dir(path);
	},

	// Create a file, given a path and string
	createFile: async (path, str) => {
		var buffer = new Uint8Array(str.length);
		buffer.set(strToChars(str));

		await get(cli).emulator.create_file(path, buffer);
	}
});
