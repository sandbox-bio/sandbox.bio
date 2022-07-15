import { merge } from "lodash";
import localforage from "localforage";
import { get, writable, derived } from "svelte/store";
import { getLocalForageKey } from "stores/config";

import awk_data from "sandboxes/orders.txt";


// =============================================================================
// Constants
// =============================================================================

// Tools to load in playground
export const TOOLS = [
	{ name: "jq", aioliConfig: { tool: "jq", version: "1.6" }},
	{ name: "awk", aioliConfig: { tool: "gawk", version: "5.1.0", reinit: true }},
	{ name: "grep", aioliConfig: { tool: "grep", version: "3.7", reinit: true }},
	{ name: "sed", aioliConfig: { tool: "sed", version: "4.8", reinit: true }}
];

// Supported flags
export const FLAG_SETTING = "setting";    // Can only use flag once (dropdown)
export const FLAG_BOOLEAN = "boolean";    // Can only use flag once; toggle between true/false (checkbox)
export const FLAG_PARAM   = "parameter";  // Can use flag multiple times (button to add new parameter)
export const FLAGS = {
	awk: [
		{
			name: "Delimiter",
			flag: "-F",
			type: FLAG_SETTING,
			values: [
				{ name: "Tab", value: `"\\t"` },
				{ name: "Comma", value: `","` },
				{ name: "Space", value: `" "` },
			]
		},
		// {
		// 	name: "Add Variable",
		// 	flag: "-v",
		// 	type: FLAG_PARAM,
		// 	value: `myVar="Some value"`
		// }
	],
	jq: [
		{
			name: "Compact",
			flag: "-c",
			type: FLAG_BOOLEAN,
		},
		{
			name: "Sort keys",
			flag: "-S",
			type: FLAG_BOOLEAN,
		},
		{
			name: "Slurp (read input into array)",
			flag: "-s",
			type: FLAG_BOOLEAN,
		}
	]
};

// Examples
export const EXAMPLES = {
	awk: [
		{
			name: "Output nth column",
			command: `{ print $3 }`,
			input: awk_data,
			flags: `-F "\\t"`
		},
		{
			name: "Filter first, then output nth column",
			command: `/Burrito/ { print $3 }`,
			input: awk_data,
			flags: `-F "\\t"`
		},
		{
			name: "Pass variables into <code>awk</code> from the outside",
			command: `# Use NR > 1 to skip the header line
NR > 1 {
	print substr($5, 2) * tax
}`,
			input: awk_data,
			flags: `-F "\\t" -v tax=0.15`
		},
	]
}

// Store defaults
const DEFAULT = {
	settings: {
		interactive: true,
	},
	data: {
		awk: {
			flags: "-F \\t",
			command: `/Burrito/ { print $3 }`,
			input: awk_data
		},
		jq: {
			flags: "-S",
			command: ".",
			input: `{"lastname": "Aboukhalil", "firstname":"Robert"}`
		}
	},
};


// =============================================================================
// Stores
// =============================================================================

// Create sandbox store
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

//
export const tool = writable(null);

// Create derived store for data
export const data = derived(
	sandbox,
	$sandbox => $sandbox?.data[get(tool)?.name]
);
