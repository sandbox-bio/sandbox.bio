// App-wide status. Keep this separate from config.js to avoid circular dependencies.
import { writable } from "svelte/store";

export const status = writable({
	app: null,
	terminal: null
});
