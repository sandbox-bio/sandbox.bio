import { get, readable, writable } from "svelte/store";
import { createClient } from "@supabase/supabase-js";
import { env } from "$env/dynamic/public";

// -----------------------------------------------------------------------------
// Variables that we'll export (_var will be exported as $var)
// -----------------------------------------------------------------------------

const _supabase = createClient(env.PUBLIC_SUPABASE_URL, env.PUBLIC_SUPABASE_API_KEY);
const _user = _supabase.auth.user();

// -----------------------------------------------------------------------------
// Exports
// -----------------------------------------------------------------------------

// User state
export const supabase = readable(_supabase);
export const user = writable(_user);
export const progress = writable({});

// // -----------------------------------------------------------------------------
// // On change
// // -----------------------------------------------------------------------------

// // Fetch information from localForage / DB
// export async function envInit() {
// 	// Refresh user. This fixes the issue where you go to sandbox.bio and it shows
// 	// you as logged out, but when you refresh it shows you as logged in.
// 	// To test:
// 	// 	  data = JSON.parse(localStorage.getItem("supabase.auth.token"))
// 	// 	  data.currentSession.expires_at = data.currentSession.expires_at - 1e6
// 	// 	  data.expiresAt = data.expiresAt - 1e6
// 	// 	  localStorage.setItem("supabase.auth.token", JSON.stringify(data))
// 	if (_supabase.auth.user() === null)
// 		_supabase.auth.onAuthStateChange((event, session) => {
// 			if (event === "SIGNED_IN") window.location.reload();
// 		});

// 	// Initialize
// 	status.set({ ...status, app: null });
// 	console.log("envInit()", get(user)?.email);
// 	let dataEnv = {},
// 		dataProgress = {};

// 	// User is logged out ==> use env vars from localForage
// 	if (get(user) === null) {
// 		dataEnv = await localforage.getItem(getLocalForageKey());

// 		// User is logged in ==> use env vars from DB
// 	} else {
// 		const data = (await _supabase.from(`state${tableAppend}`).select()).data; // always returns array, even if empty
// 		dataEnv = data[0]?.env;
// 		dataProgress = data[0]?.progress;
// 	}

// 	// Make sure default env vars are all defined
// 	dataEnv = dataEnv || {};
// 	dataProgress = dataProgress || {};
// 	for (let v in _config.env) if (!dataEnv[v]) dataEnv[v] = _config.env[v];

// 	// Update env variable but don't update DB because that's where we got the data from!
// 	await env.set(dataEnv);
// 	await progress.set(dataProgress);
// 	status.set({ ...status, app: true });
// }

// // When user logs in (object) or logs out (null)
// user.subscribe(async (userUpdated) => {
// 	if (!get(status).app) return;
// 	console.log("user.subscribe", userUpdated?.email);
// 	await envInit();
// });

// // When an environment variable is updated, update localForage+DB
// env.subscribe(async (envUpdated) => {
// 	if (!get(status).app) return;
// 	// Don't update the DB if no values have actually changed
// 	const envPrevious = await localforage.getItem(getLocalForageKey());
// 	if (JSON.stringify(envPrevious) == JSON.stringify(envUpdated)) return;
// 	console.log("env.subscribe", envUpdated);

// 	// Update localForage for both guest/logged in users
// 	await localforage.setItem(getLocalForageKey(), envUpdated);
// 	// And update state in DB if user is logged in
// 	if (get(user) !== null) await updateDB({ env: envUpdated });
// });

// // When tutorial progress is updated, update DB
// progress.subscribe(async (progressUpdated) => {
// 	if (!get(status).app) return;

// 	// Update DB if user is logged in
// 	if (get(user) !== null) {
// 		console.log("progress.subscribe", progressUpdated);
// 		await updateDB({ progress: progressUpdated });
// 	}
// });

// // -----------------------------------------------------------------------------
// // Utility functions
// // -----------------------------------------------------------------------------

// The key to use for storing information in localForage
export function getLocalForageKey(type = "env") {
	if (!["env", "studio", "fs", "ide", "history", "sandbox", "quiz"].includes(type)) throw `Unexpected type ${type}.`;

	let key = `${type}:${get(user)?.id || "guest"}`;
	if (["fs", "ide", "sandbox", "quiz"].includes(type)) key += ":";
	return key;
}

// // Update user state in DB
// async function updateDB(update) {
// 	console.log("updateDB", update);

// 	// Try to update
// 	const { data, error } = await _supabase
// 		.from(`state${tableAppend}`)
// 		.update(update)
// 		.match({ user_id: get(user).id });
// 	// If fails, do an insert
// 	if (!data) {
// 		update.user_id = get(user).id;
// 		await _supabase.from(`state${tableAppend}`).insert(update);
// 	}
// }
