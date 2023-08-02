import { get, writable } from "svelte/store";

// Current tutorial
export const cli = writable({
    // State
	emulator: null, // v86 emulator object once it's ready
	xterm: null, // xterm "term" object once it's ready
	addons: {
		serialize: null,
        fit: null
	}, // xterm addons once they're ready
	
    // Utilities
    exec: (cmd) => {
		const v86 = get(cli).emulator;
		const command = `${cmd}\n`;
		const chars = command.split("");
		chars.forEach((d) => {
			v86.bus.send("serial0-input", d.charCodeAt());
		});
	}
});
