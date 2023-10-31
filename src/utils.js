import { get } from "svelte/store";
import localforage from "localforage";
import { createClient } from "@supabase/supabase-js";
import { env } from "$env/dynamic/public";
import { LOGGING, LOGGING_DEBUG } from "$src/config";
import { user } from "$stores/user";

export const supabaseAnon = createClient(env.PUBLIC_SUPABASE_URL, env.PUBLIC_SUPABASE_API_KEY);

const STATE_FS = "fs";

export class LocalState {
	static async getKey(state, tutorial) {
		const prefix = get(user)?.email || "guest";
		return `${prefix}:${state}:${tutorial}`;
	}

	// File system
	static async getFS(tutorial) {
		if (!tutorial) return [];

		const key = await LocalState.getKey(STATE_FS, tutorial);
		return (await localforage.getItem(key)) || [];
	}

	static async setFS(tutorial, value) {
		if (!tutorial) throw "Stopped saving FS state because moved away from terminal.";
		log(LOGGING_DEBUG, "Saving FS state...");

		const key = await LocalState.getKey(STATE_FS, tutorial);
		return await localforage.setItem(key, value);
	}
}

// Get table name based on environment. Only use this for public.* tables
export function t(tableName) {
	if (env.PUBLIC_ENVIRONMENT !== "prd") return `${tableName}_stg`;
	return tableName;
}

export function log(level, ...message) {
	if (LOGGING >= level) {
		console.log(...message);
	}
}

export function strToChars(str) {
	const chars = str.split("");
	return chars.map((d) => d.charCodeAt(0));
}
