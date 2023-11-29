// App-wide status. Keep this separate from config.js to avoid circular dependencies.
import { writable } from "svelte/store";

export const status = writable({
	igv: null // Used to update IGV settings from an IGV tutorial
});
