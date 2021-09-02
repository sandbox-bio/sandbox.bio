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
let dbDontUpdate = false;  // is set to true once the env info was fetched from the DB
let appIsReady = false;    // is set to true once we get user login info on page load

// Fetch information from localForage / DB
export async function envInit()
{
	console.log("envInit()", _user.email);
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

	// Update env variable but don't update DB because that's where we got the data from!
	dbDontUpdate = true;
	await env.set(dataEnv);
}

// This is triggered when the user logs in / logs out (appIsReady means we don't trigger at page load)
user.subscribe(async userUpdated => {
	// Note that userUpdated == null is a valid state of being logged out
	if(!appIsReady) {
		appIsReady = true;
		return;
	}
	console.log("user.subscribe", userUpdated.email);
	await envInit();
});

// When an environment variable is updated, update localForage and DB
env.subscribe(async envUpdated => {
	if(!envUpdated || JSON.stringify(envUpdated) === "{}")
		return;
	console.log("env.subscribe", envUpdated);

	// Update localForage for both guest/logged in users
	await localforage.setItem(LOCAL_STORAGE_ENV(), envUpdated);
	// And update state in DB if user is logged in
	if(_user !== null)
		await updateDB({ env: envUpdated });
});

// Utility function to update user state
async function updateDB(update) {
	if(dbDontUpdate) {
		dbDontUpdate = false;
		return;
	}
	console.log("updateDB", update);

	// Try to update
	const { data, error } = await _supabase.from("state").update(update).match({ user_id: _user.id });
	// If fails, do an insert
	if(!data) {
		update.user_id = _user.id;
		await _supabase.from("state").insert(update);
	}
}
