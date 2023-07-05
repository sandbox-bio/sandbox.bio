// App-wide status. Keep this separate from config.js to avoid circular dependencies.
import { writable } from "svelte/store";

export const status = writable({
	app: null, // Set to false if app is not yet initialized; set to true once envInit() done running
	terminal: null, // Set to "execDone" whenever a command finishes running
	igv: null // Used to update IGV settings from an IGV tutorial
});
