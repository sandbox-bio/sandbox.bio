// Manage logging

import { createClient } from "@supabase/supabase-js";

// Config
const PATHS_IGNORE = ["/favicon.ico"];

// Supabase Client (variables stored as secrets in Cloudflare Worker)
const supabase = createClient(SUPABASE_URL, SUPABASE_API_KEY);

// Log an event and response
export async function logEvent(event, response) {
	const url = new URL(event.request.url);
	const pathname = url.pathname;

	// Ignore certain paths
	if(PATHS_IGNORE.includes(pathname))
		return;

	// Store path and status
	await supabase.from("logs").insert([{
		pathname,
		status: (await response).status
	}]);
}
