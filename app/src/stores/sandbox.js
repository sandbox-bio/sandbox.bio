import { writable } from "svelte/store";
import localforage from "localforage";
import { getLocalForageKey } from "stores/config";
import { merge } from "lodash";

// Store defaults
const DEFAULT = {
	settings: {
		interactive: true,
	}
};

// Create store
const { set, subscribe } = writable(DEFAULT);
export const sandbox = {
	set: value => {
		// This makes sure we always have defaults if a key was not previously defined in localforage
		// const data = Object.assign({}, DEFAULT, value);
		const data = merge({}, DEFAULT, value);

		// Save in localforage
		localforage.setItem(getLocalForageKey("sandbox"), data);

		// Save in memory
		return set(data);
	},
	init: async () => {
		set(await localforage.getItem(getLocalForageKey("sandbox")) || DEFAULT);
	},
	subscribe,
};
