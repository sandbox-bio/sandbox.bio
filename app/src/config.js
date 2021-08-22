import { readable, writable } from "svelte/store";

const URL_API = `https://${window.location.hostname != "localhost" ? window.location.hostname : "dev.sandbox.bio"}/api/v1`;

export const config = readable({
	api: URL_API,
	hostname: "sandbox"
});

export const vars = writable({
	PS1: "\\u@\\h$",
	USER: "guest"
});
