import localforage from "localforage";
import { LOGGING, LOGGING_DEBUG } from "$src/config";

export class LocalState {
	static getKeyFS(tutorial) {
		// FIXME: add user info
		return `guest:fs:${tutorial}`;
	}

	static async getFS(tutorial) {
		if (!tutorial) return [];

		const key = LocalState.getKeyFS(tutorial);
		return (await localforage.getItem(key)) || [];
	}

	static async setFS(tutorial, value) {
		if (!tutorial) throw "Warning: Stop saving FS state (moved away from terminal).";
		log(LOGGING_DEBUG, "Saving FS state...");

		const key = LocalState.getKeyFS(tutorial);
		return await localforage.setItem(key, value);
	}
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
