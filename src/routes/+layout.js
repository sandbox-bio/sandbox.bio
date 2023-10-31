import { supabaseAnon, t } from "$src/utils.js";

export const ssr = false;

// Fetch data before page loads
export const load = async () => {
	// Get user session
	const {
		data: { session }
	} = await supabaseAnon.auth.getSession();
	const user = session?.user || {};

	// Fetch tutorial progress
	let progress = {};
	if (user?.id) {
		const { data, error } = await supabaseAnon.from(t("state")).select().single();
		if (error) {
			console.error(error);
		} else {
			progress = data?.progress || {};
		}
	}

	return { user, progress };
};
