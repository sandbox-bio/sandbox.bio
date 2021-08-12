// API endpoint for tracking analytics around which steps people are currently on
addEventListener("fetch", (event) => {
	event.respondWith(handleRequest(event.request));
});

async function handleRequest(request) {
	// Ping requests must include "from" and "to"
    const url = new URL(request.url);
	const pingFrom = url.searchParams.get("from");
	const pingTo = url.searchParams.get("to");
	if(!pingFrom || !pingTo || request.method != "POST")
		return new Response("", { status: 404 });

	// Store in KV, where key = <timestamp>:<uuid>
	// For privacy, don't store IP addresses, only hashes of them
	const ip = await hash(request.headers.get("CF-Connecting-IP"));
	const key = `${new Date().getTime()}:${uuidv4()}`;
	await PING.put(key, "", { metadata: [ip, pingFrom, pingTo] });

	// Enable CORS
	return new Response("PONG");
}

// Generate a UUID (source: https://stackoverflow.com/a/2117523)
function uuidv4() {
	return ([1e7]+-1e3+-4e3+-8e3+-1e11).replace(/[018]/g, c =>
		(c ^ crypto.getRandomValues(new Uint8Array(1))[0] & 15 >> c / 4).toString(16)
	);
}

// MD5 hash (source: https://stackoverflow.com/a/64795218)
async function hash(message) {
	const msgUint8 = new TextEncoder().encode(message) // encode as (utf-8) Uint8Array
	const hashBuffer = await crypto.subtle.digest("MD5", msgUint8) // hash the message
	const hashArray = Array.from(new Uint8Array(hashBuffer)) // convert buffer to byte array
	const hashHex = hashArray.map(b => b.toString(16).padStart(2, '0')).join('') // convert bytes to hex string  

	return hashHex;
}
