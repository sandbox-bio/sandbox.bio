import localforage from "localforage";
import { get, readable, writable } from "svelte/store";
import { createClient } from "@supabase/supabase-js";

// -----------------------------------------------------------------------------
// Config
// -----------------------------------------------------------------------------

const hostname = window.location.hostname == "localhost" ? "dev.sandbox.bio" : window.location.hostname;
const urlsSupabase = {
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
// Variables that we'll export (_var will be exported as $var)
// -----------------------------------------------------------------------------

const _supabase = createClient(urlsSupabase[hostname].url, urlsSupabase[hostname].publicKey);
const _user = _supabase.auth.user();

const _config = {
	// API
	api: `https://${hostname}/api/v1`,
	// Default environment information
	hostname: "sandbox",
	playground: `\u001b[0;37m# This playground is for open-ended exploration.\n# For guided tutorials, see https://sandbox.bio/tutorials\n#\n# Example:\n#   samtools view -o test.bam /samtools/examples/toy.sam\n#   samtools index test.bam\n#   ls test.bam.bai\n#   samtools idxstats test.bam  # idxstats uses the .bai file \u001b[0m\n`,
	// Default environment variables. These are auto-regenerated in Terminal.svelte:input() even if the user deletes them
	env: {
		PS1: '\\u@\\h$ ',
		USER: "guest",
		HOME: "/shared/data"
	}
};

// -----------------------------------------------------------------------------
// Exports
// -----------------------------------------------------------------------------

// Login/Signup modal open or not
export const loginModal = writable(false);

// App settings (read-only)
export const config = readable(_config);

// User state
export const supabase = readable(_supabase);
export const user = writable(_user);
export const env = writable({});
export const progress = writable({});

// -----------------------------------------------------------------------------
// On change
// -----------------------------------------------------------------------------

const LOCAL_STORAGE_ENV = () => `env:${_user?.id || "guest"}`;

// Fetch information from localForage / DB
export async function initEnv()
{
	console.log("initEnv()")
	let dataEnv = {};

	// User is logged out ==> use env vars from localForage
	if(_user === null) {
		dataEnv = await localforage.getItem(LOCAL_STORAGE_ENV());

	// User is logged in ==> use env vars from DB
	} else {
		const data = (await _supabase.from("state").select()).data;  // always returns array, even if empty
		dataEnv = data[0]?.env;
	}

	// Make sure default env vars are all defined
	dataEnv = dataEnv || {};
	for(let v in _config.env)
		if(!dataEnv[v])
			dataEnv[v] = _config.env[v];

	env.set(dataEnv);
}

// This is triggered when the page loads or when the user logs in / logs out
user.subscribe(async userUpdated => {
	console.log("user.subscribe", userUpdated);
	// let dataEnv = {}, dataProgress = {};
	// progress.set(dataProgress);
});

// When an environment variable is updated, update localForage and DB
env.subscribe(async envUpdated => {
	if(JSON.stringify(envUpdated) === "{}")
		return;
	console.log("env.subscribe", envUpdated);

	// // Update localForage for both guest/logged in users
	// await localforage.setItem(LOCAL_STORAGE_ENV(), envUpdated);

	// // And update state in DB if user is logged in
	// if(_user !== null)
	// 	await updateState({ env: envUpdated });
});

