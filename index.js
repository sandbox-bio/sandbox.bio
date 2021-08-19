// Entrypoint for sandbox.bio; based on the Cloudflare Workers Site template.

import { Router } from "itty-router";
import { getAssetFromKV } from "@cloudflare/kv-asset-handler";


// =============================================================================
// Config
// =============================================================================

const routerAPI = Router({ base: "/api/v1" });
const CACHE_CONTROL = {
	browserTTL: 604800,  // 1 week (default: null)
	edgeTTL: 172800,     // 2 days (default: 2 days)
	bypassCache: false   // Do not bypass Cloudflare's cache (default: false)
};


// =============================================================================
// Routes
// =============================================================================

routerAPI.get('/todos', () => new Response('Todos Index!'));
routerAPI.get('/todos/:id', ({ params }) => new Response(`Todo #${params.id}`));
routerAPI.post('/todos', async request => {
	const content = await request.json();
	return new Response('Creating Todo: ' + JSON.stringify(content));
})
routerAPI.all('*', () => new Response('Not Found.', { status: 404 }));


// =============================================================================
// Process requests
// =============================================================================

addEventListener("fetch", event => {
	let response = {};

	// Process user request
	// TODO: Add logging for successes and failures
	try {
		const url = new URL(event.request.url);
		if(url.pathname.startsWith("/api/"))
			response = routerAPI.handle(event.request);
		else
			response = handleStaticPage(event);
	} catch (e) {
		response = new Response("Internal Error", { status: 500 });
	}

	// Return result
	event.respondWith(response);
})

async function handleStaticPage(event) {
	// Retrieve KV value given a path, or show 404
	try {
		return await getAssetFromKV(event, {
			mapRequestToAsset: handlePrefix,
			cacheControl: CACHE_CONTROL
		});
	} catch (e) {
		return new Response("404 Not Found", { status: 404 });
	}
}

function handlePrefix(request) {
	// If path looks like a directory, append index.html
	const parsedUrl = new URL(request.url);
	let pathname = parsedUrl.pathname;
	if (pathname.endsWith("/"))
		pathname = pathname.concat("index.html");

	parsedUrl.pathname = pathname;
	return new Request(parsedUrl.toString(), request);
}
