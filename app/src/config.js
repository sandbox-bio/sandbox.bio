import { readable, writable } from "svelte/store";

// App settings (read-only)
export const config = readable({
	// API
	api: `https://${window.location.hostname != "localhost" ? window.location.hostname : "dev.sandbox.bio"}/api/v1`,
	// Default environment information
	hostname: "sandbox",
	env: {
		PS1: '\\u@\\h$ ',
		USER: "guest"
	}
});

// User-defined variables
export const vars = writable({});
