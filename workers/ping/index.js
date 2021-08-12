// API endpoint for tracking analytics around which steps people are currently on
addEventListener("fetch", (event) => {
	event.respondWith(handleRequest(event.request));
});

async function handleRequest(request) {
	// Ping requests must include "from" and "to"
    const url = new URL(request.url);
	const pingFrom = url.searchParams.get("from");
	const pingTo = url.searchParams.get("to");
	if(!pingFrom || !pingTo)
		return new Response("", { status: 404 });

	// Store in KV, where key = <timestamp>:<uuid>
	const key = `${new Date().getTime()}:${uuidv4()}`;
	await PING.put(key, "", { metadata: [pingFrom, pingTo] });

	// Enable CORS
	return new Response("PONG");
}

// Generate a UUID (source: https://stackoverflow.com/a/2117523)
function uuidv4() {
	return ([1e7]+-1e3+-4e3+-8e3+-1e11).replace(/[018]/g, c =>
		(c ^ crypto.getRandomValues(new Uint8Array(1))[0] & 15 >> c / 4).toString(16)
	);
}
