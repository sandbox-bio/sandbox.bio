export async function POST({ request }) {
	const data = await request.json();

	return new Response("pong");
}
