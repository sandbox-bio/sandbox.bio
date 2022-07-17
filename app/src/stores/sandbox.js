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
			name: "Input Delimiter",
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
			name: "Output columns",
			input: awk_data,
			flags: `-F "\\t"`,
			command: `{
	print $3, $5
}`
		},
		{
			name: "Filter and output columns",
			input: awk_data,
			flags: `-F "\\t"`,
			command: `/Burrito/ {
	print $3, $5
}`
		},
		{
			name: "Sum over a column",
			input: awk_data,
			flags: `-F "\\t"`,
			command: `BEGIN {
	sum = 0     # The BEGIN block is optional
} {
	sum += $2
} END {
	print(sum)  # The END block is executed at the very end
}`
		},
		{
			name: "Calculate a new column",
			input: awk_data,
			flags: `-F "\\t"`,
			command: `BEGIN {
	tax = 0.0725
	OFS = "\\t"  # Use tab as output separator
} {
	if (NR == 1) {
		# Rename header
		$6 = "item_price_plus_tax"
    } else {
		# Remove dollar sign from price in 5th column
		price = substr($5, 2)

		# Add 6th column with price including tax
		$6 = sprintf("$%.2f", price * (1 + tax))
	}

	# Output all columns
	print
}`
		},
		{
			name: "Inject variables",
			input: awk_data,
			flags: `-F "\\t" -v food=Chicken`,
			command: `{
	# Only print lines that match the food variable defined in the flags
	# Try changing the variable to "Tacos"
	if($3 ~ food)
		print
}`,
		},
		{
			name: "Dictionaries and loops",
			input: awk_data,
			flags: `-F "\\t"`,
			command: `# Skip first line (header)
NR > 1 {
	itemCount = $2
	itemName = $3

	# Track how often burritos were ordered
	# Notice that counts becomes a dictionary with
	# default values of zero automatically!
	if(itemName ~ /Burrito/)
		counts[itemName] += itemCount
} END {
	# Print burrito counts, split by filling
	for(item in counts)
		print(item, counts[item])
}`,
		},
		{
			name: "Cumulative sums on groups of rows",
			input: awk_data,
			flags: `-F "\\t"`,
			command: `BEGIN {
	currentOrderId = -1
	totalPrice = 0

# Print header
} NR == 1 {
	print "order_id", "order_price"

# Add up prices of all items within the same order
# (assumes input is sorted by the order_id column)
} NR > 1 {
	# If we come across a new order ID, output the
	# total price of the previous order.
	if($1 != currentOrderId) {
		if(currentOrderId != -1)
			print currentOrderId, totalPrice

		# Reset price
		currentOrderId = $1
		totalPrice = 0
    }

	totalPrice += substr($5, 2)

# Make sure you don't forget the last order in the file!
} END {
	print currentOrderId, totalPrice
}`
		},
		{
			name: "Functions",
			input: awk_data,
			flags: `-F "\\t"`,
			command: `# Sanitize data so it's friendlier for analysis in other tools
function sanitizeStr(str) {
	# Replace spaces with underscores (\\s = whitespace)
	gsub(/\\s/, "_", str)

	# Remove square brackets
	gsub(/\\[|\\]/, "", str)

	return str

} BEGIN {
	OFS = "\\t"

} {
	# Overwrite 3rd/4th columns (doesn't affect source files!)
	$3 = sanitizeStr($3)
	$4 = sanitizeStr($4)
	print
}`
		}
	],
	jq: []
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
