import { env } from "$env/dynamic/private";
import { json } from "@sveltejs/kit";

export async function POST({ request }) {
	const data = await request.json();
	const payload = {
		method: "POST",
		headers: {
			"Content-Type": "application/json",
			Accept: "application/json",
			"X-Kit-Api-Key": env.KIT_API_KEY
		},
		body: JSON.stringify({
			email_address: data.email
		})
	};

	// Add email to subscribes
	await fetch(`https://api.kit.com/v4/subscribers`, payload).then((d) => d.json());

	// Tag email with sandbox.bio
	const response = await fetch(`https://api.kit.com/v4/tags/${env.KIT_TAG_ID}/subscribers`, payload).then((d) => d.json());

	return json({
		error: response?.errors?.join("\n")
	});
}
