// FIXME: use LocalState from $stores/user instead
import { get } from "svelte/store";

// The key to use for storing information in localForage
export function getLocalForageKey(type = "env") {
	if (!["env", "studio", "fs", "ide", "history", "sandbox", "quiz"].includes(type)) throw `Unexpected type ${type}.`;

	let key = `${type}:${get(user)?.id || "guest"}`;
	if (["fs", "ide", "sandbox", "quiz"].includes(type)) key += ":";
	return key;
}
