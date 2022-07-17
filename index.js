// Entrypoint for sandbox.bio; based on the Cloudflare Workers Site template.
// TODO: Add logging for successes and failures in API calls

import { logEvent } from "./api/logs";
import { routerAPI } from "./api/index";
import { getAssetFromKV } from "@cloudflare/kv-asset-handler";

// =============================================================================
// Config
// =============================================================================

// Paths to reroute to index.html so Svelte can do the routing
const PATHS_ROUTE_TO_MAIN = ["/tutorials", "/playground", "/rosalind", "/sandboxes"];

// Cache settings for static content (API is not cached)
const CACHE_CONTROL = {
	edgeTTL    : 172800,  // 2 days (default: 2 days)
	browserTTL : 604800,  // 1 week (default: null)
	bypassCache:  false,  // Do not bypass Cloudflare's cache (default: false)
};

// =============================================================================
// Process requests
// =============================================================================

// Process all incoming requests to [stg.]sandbox.bio/*
addEventListener("fetch", event => {
	let response = {};

	// Process user request
	try {
		const url = new URL(event.request.url);
		// Use separate router for API
		if(url.pathname.startsWith("/api/"))
			response = routerAPI.handle(event.request);
		// But use Cloudflare Worker Sites logic for everything else
		else
			response = handleStaticPage(event);
	} catch (e) {
		response = new Response("Internal Error", { status: 500 });
	}

	// Log event *after* we return the response to the user
	event.waitUntil(logEvent(event, response));

	// Return result
	event.respondWith(response);
})

// Retrieve KV value given a path, or show 404
async function handleStaticPage(event) {
	try {
		return await getAssetFromKV(event, {
			mapRequestToAsset: handlePrefix,
			cacheControl: CACHE_CONTROL
		});
	} catch (e) {
		return new Response("404 Not Found", { status: 404 });
	}
}

// Custom logic for handling directories and routes
function handlePrefix(request) {
	const parsedUrl = new URL(request.url);
	let pathname = parsedUrl.pathname;

	// Route app paths to index.html so that Svelte can do the routing
	for(let path of PATHS_ROUTE_TO_MAIN)
		if(pathname == path || pathname.startsWith(`${path}/`))
			pathname = "index.html";

	// If path looks like a directory, append index.html
	else if (pathname.endsWith("/"))
		pathname = pathname.concat("index.html");

	parsedUrl.pathname = pathname;
	return new Request(parsedUrl.toString(), request);
}
