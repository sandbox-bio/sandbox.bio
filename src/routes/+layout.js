import { supabaseAnon, t } from "$src/utils.js";
import { user } from "$stores/user";

export const ssr = false;

// Fetch data before page loads
export const load = async () => {
	// This fixes the issue where you get redirected to a blank page after logging in
	supabaseAnon.auth.onAuthStateChange(async (event) => {
		if (event === "SIGNED_IN") {
			window.location.reload();
			return {};
		}
	});

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
