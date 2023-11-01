import { supabaseAnon, t } from "$src/utils.js";
import { user } from "$stores/user";

// Fetch data before page loads
export const load = async () => {
	// Get user session
	const {
		data: { session }
	} = await supabaseAnon.auth.getSession();
	const userInfo = session?.user || {};

	// Fetch tutorial progress
	let progress = {};
	if (userInfo?.id) {
		const { data, error } = await supabaseAnon.from(t("state")).select().single();
		if (error) {
			console.error(error);
		} else {
			progress = data?.progress || {};
		}
	}
	user.set(userInfo);

	return { progress };
};
