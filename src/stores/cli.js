import { get, writable } from "svelte/store";

// Current tutorial
export const cli = writable({
	v86: null, // V86 object once it's ready
	exec: (cmd) => {
		const v86 = get(cli).v86;
		const command = `${cmd}\n`;
		const chars = command.split("");
		chars.forEach((d) => {
			v86.bus.send("serial0-input", d.charCodeAt());
		});
	}
});
