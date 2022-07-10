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

// Tools to load in playground
export const TOOLS = [
	{ name: "jq", aioliConfig: { tool: "jq", version: "1.6" }},
	{ name: "awk", aioliConfig: { tool: "gawk", version: "5.1.0", reinit: true }},
	{ name: "grep", aioliConfig: { tool: "grep", version: "3.7", reinit: true }},
	{ name: "sed", aioliConfig: { tool: "sed", version: "4.8", reinit: true }}
];

// Supported flags
export const FLAGS = {
	awk: [
		{
			name: "Set delimiter",
			options: [
				{ name: "Tabs", flag: "-F", value: "\\t" },
				{ name: "Commas", flag: "-F", value: "," },
				{ name: "Spaces", flag: "-F", value: `" "` }
			]
		},
		{
			name: "Variables",
			options: [
				{ name: "Add new variable", flag: "-v", value: `myVar="Hello World"`, multiple: true }
			]
		}
	],
	jq: [
		{
			name: "Output",
			options: [
				{ name: "Compact view", flag: "-c", boolean: true },
				{ name: "Sorted keys", flag: "-S", boolean: true },
				{ name: "Tab indentation", flag: "--tab", boolean: true }
			]
		},
		{
			name: "Input",
			options: [
				{ name: "Read input into array (slurp)", flag: "-s", boolean: true }
			]
		},
		{
			name: "Variables",
			options: [
				{ name: "Add new variable", flag: "--arg", value: `myVar "Hello World"`, multiple: true }
			]
		}
	]
};
