// // The biowasm CDN is hosted using Cloudflare Worker Sites where CDN files
// // are stored in a key-value store (Cloudflare Worker KV). For example:
// //     key=samtools/1.10/samtools.wasm --> value=<samtools.wasm contents>
// // 
// // A Cloudflare Worker (a serverless function) is the entry point for
// // retrieving those files from the key-value store. This index.js file
// // defines the code for that entry point. The code is mostly using the
// // Cloudflare Workers Site template, but modified to enable CORS (so that
// // biowasm modules can be loaded from non-biowasm.com domains!) and to
// // log basic stats about the number of downloads per module.

import { getAssetFromKV } from "@cloudflare/kv-asset-handler";

addEventListener("fetch", event => {
	let response = {};

	// Process user request
	try {
		// TODO: Do routing with /app/playground
		// TODO: Do routing with /app/tutorials/<id>

		// TODO: Differntiate app and home page
		// if(path.startsWith("/app"))
		// if(path == "/")
		response = handleStaticPage(event);

		// if(path.startsWith("/api/v1/"))
		// TODO: API
	} catch (e) {
		response = new Response("Internal Error", { status: 500 });
	}

	// Log basic stats
	// TODO:

	// Return result
	event.respondWith(response);
})

async function handleStaticPage(event) {
	// Retrieve KV value given a path
	try {
		let response = await getAssetFromKV(event, {
			mapRequestToAsset: handlePrefix(),
			cacheControl: {
				browserTTL: 604800,  // 1 week (default: null)
				edgeTTL: 172800,     // 2 days (default: 2 days)
				bypassCache: false   // Do not bypass Cloudflare's cache (default: false)
			}
		});

		// // Enable CORS
		// response.headers.set("Access-Control-Allow-Origin", "*");
		return response;
	// On error, show a 404 error
	} catch (e) {
		return new Response("404 Not Found", { status: 404 });
	}
}

function handlePrefix() {
	return request => {
		// If path looks like a directory append index.html
		// e.g. If path is /about/ -> /about/index.html
		const parsedUrl = new URL(request.url);
		let pathname = parsedUrl.pathname;
		if (pathname.endsWith("/"))
			pathname = pathname.concat("index.html");

		parsedUrl.pathname = pathname;
		return new Request(parsedUrl.toString(), request);
	}
}

