import { env } from "$env/dynamic/private";
import { json } from "@sveltejs/kit";

export async function POST({ request }) {
	const data = await request.json();
	const response = await fetch(`https://api.convertkit.com/v3/tags/${env.CONVERTKIT_TAG_ID}/subscribe`, {
		method: "POST",
		headers: { "Content-Type": "application/json; charset=utf-8" },
		body: JSON.stringify({
			api_secret: env.CONVERTKIT_API_KEY,
			email: data.email
		})
	}).then((d) => d.json());

	return json({
		error: response.message
	});
}
