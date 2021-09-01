import localforage from "localforage";
import { get, readable, writable } from "svelte/store";
import { createClient } from "@supabase/supabase-js";

// -----------------------------------------------------------------------------
// Config
// -----------------------------------------------------------------------------

const hostname = window.location.hostname == "localhost" ? "dev.sandbox.bio" : window.location.hostname;
const urlAPI = `https://${hostname}/api/v1`;
const urlSupabase = {
	"dev.sandbox.bio": {
		url: "https://bqjvxpdzkembvixymfae.supabase.co",
		publicKey: "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJyb2xlIjoiYW5vbiIsImlhdCI6MTYyOTMyODIwNCwiZXhwIjoxOTQ0OTA0MjA0fQ.7DzKM4bOGK1t-pPkfSe-2ALxcW5xWwcsaZfbCMWDBbY"
	},
	"stg.sandbox.bio": {
		url: "",
		publicKey: ""
	},
	"prd.sandbox.bio": {
		url: "",
		publicKey: ""
	}
};

// -----------------------------------------------------------------------------
// Exports
// -----------------------------------------------------------------------------

// User-defined variables
export const env = writable({});

// Supabase client
export const supabase = readable(createClient(urlSupabase[hostname].url, urlSupabase[hostname].publicKey));
export const user = writable(get(supabase).auth.user());

// App settings (read-only)
export const config = readable({
	// API
	api: urlAPI,
	// Default environment information
	hostname: "sandbox",
	// Default environment variables. These are auto-regenerated in Terminal.svelte:input() even if the user deletes them
	env: {
		PS1: '\\u@\\h$ ',
		USER: "guest",
		HOME: "/shared/data"
	},
	playground: `\u001b[0;37m# This playground is for open-ended exploration.
# For guided tutorials, see https://sandbox.bio/tutorials
#
# Example:
#   samtools view -o test.bam /samtools/examples/toy.sam
#   samtools index test.bam
#   ls test.bam.bai
#   samtools idxstats test.bam  # idxstats uses the .bai file \u001b[0m
`
});

// -----------------------------------------------------------------------------
// On change
// -----------------------------------------------------------------------------

const LOCAL_STORAGE_ENV = () => `env:${get(user)?.id || "guest"}`;

// This is triggered when the page loads or when the user logs in / logs out
user.subscribe(async userUpdated => {
	let dataEnv = {};

	// User is logged out ==> use env vars from local storage
	if(userUpdated === null) {
		dataEnv = await localforage.getItem(LOCAL_STORAGE_ENV());

	// User is logged in ==> use env vars from DB
	} else {
		dataEnv = (await get(supabase).from("state").select()).data;
		if(dataEnv?.length == 0)
			dataEnv = null;
		else
			dataEnv = dataEnv[0]?.env;
	}

	// If no data, initialize with default env vars
	if(dataEnv === null)
		dataEnv = get(config).env;
	env.set(dataEnv);
});

env.subscribe(async envUpdated => {
	if(JSON.stringify(envUpdated) === "{}")
		return;

	// Update localStorage for both guest/logged in users
	await localforage.setItem(LOCAL_STORAGE_ENV(), envUpdated);

	// And update state in DB if user is logged in
	if(get(user) !== null) {
		const { data, error } = await get(supabase)
			.from("state")
			.update({ env: envUpdated })
			.match({ user_id: get(user).id });
		if(!data) {
			const { data, error } = await get(supabase)
				.from("state")
				.insert({ env: envUpdated, user_id: get(user).id });
		}
	}
});
