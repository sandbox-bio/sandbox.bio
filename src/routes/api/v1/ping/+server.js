import { supabase, hash, t } from "$routes/api/utils";

export async function POST({ request }) {
	const data = await request.json();

    await supabase.from(t("pings")).insert([{
		ip: await hash(request.headers.get("CF-Connecting-IP")) || "local",
		tutorial: data.tutorial,
		step_from: data.from ? +data.from : null,
		step_to: data.to ? +data.to : null,
		playground: data.playground,
		example: data.example
	}]);

	return new Response("pong");
}
