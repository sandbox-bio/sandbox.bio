// Manage logging

import { supabase, t } from "./utils";

// =============================================================================
// Config
// =============================================================================

const PATHS_INCLUDE = ["/", "/rosalind", "/rosalind/", "/playground", "/playground/", "/playgrounds", "/playgrounds/", "/tutorials", "/tutorials/"];


// =============================================================================
// Utility functions
// =============================================================================

// Log an event and response
export async function logEvent(event, response) {
	const url = new URL(event.request.url);
	const pathname = url.pathname;
	const search = url.search;

	// Ignore certain paths
	if(!PATHS_INCLUDE.includes(pathname))
		return;

	// Store path and status
	await supabase.from(t("logs")).insert([{
		pathname,
		search,
		status: (await response).status
	}]);
}
