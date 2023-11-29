import { get, writable } from "svelte/store";
import { browser } from "$app/environment";
import { supabaseAnon, t } from "$src/utils";
import { user } from "$stores/user";

export const progress = writable({});

// When tutorial progress is updated, update DB
progress.subscribe(async (newProgress) => {
	const userInfo = get(user);
	if (!browser || !userInfo.id) return;
	const query = supabaseAnon.from(t("state")).update({ progress: newProgress }).match({ user_id: userInfo.id });

	const { error } = await query;
	// If fails, do an insert
	if (error) {
		console.error(error);
		const { error: error2 } = await supabaseAnon.from(t("state")).insert({ user_id: userInfo.id, progress: newProgress });
		console.error(error2);
	}
});
