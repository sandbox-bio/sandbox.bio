import { get, writable } from "svelte/store";
import { browser } from "$app/environment";

import { user } from "$stores/user";

export const progress = writable({});

// When tutorial progress is updated, update DB
progress.subscribe(async (newProgress) => {
	const userInfo = get(user);
	if (!browser || !userInfo.id) return;
});
