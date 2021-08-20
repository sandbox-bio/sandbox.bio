// Script I ran to migrate the ping.sandbox.bio KV data to Supabase.
// Note that the code that pushes data to the DB and deletes keys is commented out here.

addEventListener("fetch", (event) => {
	event.respondWith(
		handleRequest(event.request).catch(
			(err) => new Response(err.stack, { status: 500 })
		)
	);
});

async function handleRequest(request) {
	const url = new URL(request.url);

	// Get all keys we want to process
	let hosts = {};
	let dataCurr = {};
	while(!dataCurr.list_complete)
	{
		let dataToPush = [];
		let toDelete = [];
		dataCurr = await KV.list({ cursor: dataCurr.cursor });

		console.log(`Processing ${dataCurr.keys.length} items...`)
		dataCurr.keys.map(d => {
			const urlFrom = new URL(d.metadata[1]);
			const urlTo = new URL(d.metadata[2]);

			if(!(urlFrom.hostname in hosts))
			hosts[urlFrom.hostname] = 0;
			hosts[urlFrom.hostname]++;

			if(urlFrom.hostname != "sandbox.bio")
			{
				toDelete.push(d.name);
				dataToPush.push({
					time: new Date(+d.name.split(":")[0]),
					ip: d.metadata[0],
					tutorial: urlFrom.searchParams.get("id"),
					step_from: +urlFrom.searchParams.get("step"),
					step_to: +urlTo.searchParams.get("step"),
				});
			}
		});

		// Push data
		if(dataToPush.length > 0)
		{
			// const result = await fetch(`${SUPABASE_URL}/rest/v1/pings`, {
			// 	method: "POST",
			// 	headers: {
			// 		"apikey": SUPABASE_API_KEY,
			// 		"Content-Type": "application/json",
			// 		"Authorization": `Bearer ${SUPABASE_API_KEY}`,
			// 	},
			// 	body: JSON.stringify(dataToPush)
			// });

			if(result.status == 201)
			{
				// for(let key of toDelete)
				// 	await KV.delete(key);
			}
		}
	}

	return new Response(JSON.stringify("done"), {
		headers: { "Content-Type": "application/json" },
	});
}