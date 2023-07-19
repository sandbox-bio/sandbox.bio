import { merge } from "lodash";
import localforage from "localforage";
import { writable } from "svelte/store";
import { getLocalForageKey } from "$stores/config";

import data_text from "$components/playgrounds/orders.tsv?raw";
import data_json from "$components/playgrounds/orders.json?raw";

// =============================================================================
// Constants
// =============================================================================

// Tools to load in playground
export const TOOLS = [
	{ name: "jq", aioliConfig: { tool: "jq", version: "1.6", reinit: true } },
	{ name: "awk", aioliConfig: { tool: "gawk", version: "5.1.0", reinit: true } },
	{ name: "grep", aioliConfig: { tool: "grep", version: "3.7", reinit: true } },
	{ name: "sed", aioliConfig: { tool: "sed", version: "4.8", reinit: true } }
];

// Supported flags
export const FLAG_SETTING = "setting"; // Can only use flag once (dropdown)
export const FLAG_BOOLEAN = "boolean"; // Can only use flag once; toggle between true/false (checkbox)
export const FLAG_PARAM = "parameter"; // Can use flag multiple times (button to add new parameter)
export const FLAGS = {
	awk: [
		{
			name: "Input Delimiter",
			flag: "-F",
			type: FLAG_SETTING,
			values: [
				{ name: "Tab", value: `"\\t"` },
				{ name: "Comma", value: `","` },
				{ name: "Space", value: `" "` }
			]
		}
		// {
		// 	name: "Add Variable",
		// 	flag: "-v",
		// 	type: FLAG_PARAM,
		// 	value: `myVar="Some value"`
		// }
	],
	grep: [
		{
			name: "Invert match",
			flag: "-v",
			type: FLAG_BOOLEAN
		},
		{
			name: "Regex",
			flag: "-E",
			type: FLAG_BOOLEAN
		},
		{
			name: "Ignore case",
			flag: "-i",
			type: FLAG_BOOLEAN
		},
		{
			name: "Count",
			flag: "-c",
			type: FLAG_BOOLEAN
		},
		{
			name: "Line numbers",
			flag: "-n",
			type: FLAG_BOOLEAN
		}
	],
	sed: [
		{
			name: "Regex",
			flag: "-E",
			type: FLAG_BOOLEAN
		}
	],
	jq: [
		{
			name: "Compact",
			flag: "-c",
			type: FLAG_BOOLEAN
		},
		{
			name: "Sort keys",
			description: "Sort keys in the output",
			flag: "-S",
			type: FLAG_BOOLEAN
		},
		{
			name: "Slurp",
			description: "Read input into array",
			flag: "-s",
			type: FLAG_BOOLEAN
		},
		{
			name: "Raw input",
			flag: "-R",
			type: FLAG_BOOLEAN
		}
	]
};

// Examples
export const EXAMPLES = {
	awk: [
		{
			name: "Output columns",
			input: data_text,
			flags: `-F "\\t"`,
			command: `{
	print $3, $5
}`
		},
		{
			name: "Filter and output columns",
			input: data_text,
			flags: `-F "\\t"`,
			command: `/Burrito/ {
	print $3, $5
}`
		},
		{
			name: "Sum over a column",
			input: data_text,
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
			input: data_text,
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
			input: data_text,
			flags: `-F "\\t" -v food=Chicken`,
			command: `{
	# Only print lines that match the food variable defined in the flags
	# Try changing the variable to "Tacos"
	if($3 ~ food)
		print
}`
		},
		{
			name: "Dictionaries and loops",
			input: data_text,
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
}`
		},
		{
			name: "Cumulative sums on groups of rows",
			input: data_text,
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
			print currentOrderId, "$" totalPrice

		# Reset price
		currentOrderId = $1
		totalPrice = 0
    }

	totalPrice += substr($5, 2)

# Make sure you don't forget the last order in the file!
} END {
	print currentOrderId, "$" totalPrice
}`
		},
		{
			name: "Functions",
			input: data_text,
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
	grep: [
		{
			name: "Basic filter",
			input: data_text,
			flags: ``,
			command: `Burrito`
		},
		{
			name: "Reverse filter",
			input: data_text,
			flags: `-v`,
			command: `NULL`
		},
		{
			name: "Case insensitive filter",
			input: data_text,
			flags: `-i`,
			command: `bowl`
		},
		{
			name: "Simple regular expressions",
			input: data_text,
			flags: `-E`,
			command: `(Chicken|Carnitas) Burrito`
		},
		{
			name: "String patterns",
			input: data_text,
			flags: `-E`,
			command: `Chicken.*Fajita`
		},
		{
			name: "Lines that start with a pattern",
			input: data_text,
			flags: `-E`,
			command: `^1[0-9]`
		},
		{
			name: "Count matches",
			input: data_text,
			flags: `-c`,
			command: `Burrito`
		},
		{
			name: "Output line numbers",
			input: data_text,
			flags: `-n`,
			command: `Chips and Guacamole`
		},
		{
			name: "Show lines before/after match",
			input: data_text,
			flags: `-A 1 -B 1 --group-separator "=====" -n`,
			command: `Canned Soda`
		}
	],
	sed: [
		{
			name: "Replace first occurence on each line",
			input: data_text,
			flags: ``,
			command: `s/a/*/`
		},
		{
			name: "Replace all occurences",
			input: data_text,
			flags: ``,
			command: `s/a/*/g`
		},
		{
			name: "Replace case insensitive",
			input: data_text,
			flags: ``,
			command: `s/a/*/gi`
		},
		{
			name: "Replace with regex",
			input: data_text,
			flags: `-E`,
			command: `s/\\w/*/g`
		},
		{
			name: "Replace price pattern",
			input: data_text,
			flags: `-E`,
			command: `s/\\$[0-9]+.[0-9]+/Price Redacted/g`
		},
		{
			name: "Replace everything in brackets",
			input: data_text,
			flags: `-E`,
			command: `s/\\[.*\\]/*/g`
		},
		{
			name: "Extract subset of rows",
			input: data_text,
			flags: `-n`,
			command: `1,5p`
		},
		{
			name: "Remove header line",
			input: data_text,
			flags: `-E`,
			command: `1d`
		},
		{
			name: "Remove lines using a pattern",
			input: data_text,
			flags: `-E`,
			command: `/Burrito|Tacos/d`
		},
		{
			name: "Convert tabs to commas (TSV to CSV)",
			input: data_text,
			flags: `-E`,
			command: `s/,//g;		# Remove existing commas from description column
s/\\t/,/g;	# Convert to CSV
`
		},
		{
			name: "Multiple operations",
			input: data_text,
			flags: `-E`,
			command: `/Burrito|Tacos|Bowl/d;	# Remove entrées
1d;						# Remove header
s/\\$/USD/g;				# Replace dollar sign with USD
s/,//g					# Before converting to CSV, remove existing commas
s/\\t/,/g;				# Convert to CSV file
s/NULL//g;				# Remove NULLs
s/\\[|\\]//g				# Remove brackets
`
		}
	],
	jq: [
		{
			name: "Output formatted JSON",
			input: data_json,
			flags: ``,
			command: `.`
		},
		{
			name: "Output compact JSON",
			input: data_json,
			flags: `-c`,
			command: `.`
		},
		{
			name: "Sort keys in output",
			input: data_json,
			flags: `-S`,
			command: `.`
		},
		{
			name: "Count array elements",
			input: data_json,
			flags: ``,
			command: `. | length`
		},
		{
			name: "Get one order name",
			input: data_json,
			flags: ``,
			command: `.[5].name`
		},
		{
			name: "Get all order names",
			input: data_json,
			flags: ``,
			command: `.[] | .name`
		},
		{
			name: "Get all order names as array",
			input: data_json,
			flags: ``,
			command: `[ .[] | .name ]`
		},
		{
			name: "Get all order names as array using map()",
			input: data_json,
			flags: ``,
			command: `map(.name)`
		},
		{
			name: "Filter by value",
			input: data_json,
			flags: ``,
			command: `map(select(.name == "Chicken Bowl"))`
		},
		{
			name: "Filter by value and extract attributes",
			input: data_json,
			flags: ``,
			command: `map(select(.name == "Chicken Bowl") | .id)`
		},
		{
			name: "Filter by value in a nested array",
			input: data_json,
			flags: ``,
			command: `map(select(.description | contains(["Cheese", "Pinto Beans"]) ))`
		},
		{
			name: "Rename fields",
			input: data_json,
			flags: ``,
			command: `map({
  order_id: .id,
  item_price: .price,
  name,
  quantity,
  description,
})
`
		},
		{
			name: "Convert a TSV file to JSON",
			input: data_text,
			flags: `-R -s`,
			command: `# Ignore first and last line (header and empty line)
split("\\n")[1:-1] |

	# Split fields by tab. This generates an array,
	# accessible by .[0], .[1], etc.
	map(split("\\t")) |

	# Define a dictionary for order using array elements
	map({
		"id": .[0] | tonumber,
		"quantity": .[1] | tonumber,
		"name": .[2],
		"price": .[4] | rtrimstr(" "),
		"description": .[3] |
			split(",") |
			map(
      			ltrimstr(" ") | rtrimstr(" ") | gsub("\\\\[|]"; "")
      		)
	})
`
		}
	]
};

// Store defaults
const DEFAULT = {
	data: {},
	settings: {
		interactive: true
	}
};

TOOLS.forEach((t) => {
	DEFAULT.data[t.name] = {
		flags: EXAMPLES[t.name]?.[0]?.flags,
		command: EXAMPLES[t.name]?.[0]?.command,
		input: EXAMPLES[t.name]?.[0]?.input
	};
});

// =============================================================================
// Stores
// =============================================================================

// Create sandbox store
const { set, subscribe } = writable(DEFAULT);
export const sandbox = {
	set: (value) => {
		// This makes sure we always have defaults if a key was not previously defined in localforage
		// const data = Object.assign({}, DEFAULT, value);
		const data = merge({}, DEFAULT, value);

		// Save in localforage
		localforage.setItem(getLocalForageKey("sandbox"), data);

		// Save in memory
		return set(data);
	},
	init: async () => {
		set((await localforage.getItem(getLocalForageKey("sandbox"))) || DEFAULT);
	},
	subscribe
};
