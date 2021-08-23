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
	},
	playground: `\u001b[0;37m# This playground is for open-ended exploration.
# For guided tutorials, see https://sandbox.bio/tutorials
#
# Example:
#   samtools view -o test.bam /samtools/examples/toy.sam
#   samtools index test.bam
#   ls test.bam.bai
#   samtools idxstats test.bam  # idxstats uses the .bai file \u001b[0m`
});

// User-defined variables
export const vars = writable({});
