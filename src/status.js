// App-level status. Keep this separate from config.js to avoid circular dependencies.
import { writable } from "svelte/store";

export const status = writable({
	terminal: "ready"
});
