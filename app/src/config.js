import { readable, writable } from "svelte/store";

// App settings (read-only)
export const config = readable({
	hostname: "sandbox",
	api: `https://${window.location.hostname != "localhost" ? window.location.hostname : "dev.sandbox.bio"}/api/v1`
});

// User-defined variables
export const vars = writable({
	PS1: "\\u@\\h$ ",
	USER: "guest"
});
