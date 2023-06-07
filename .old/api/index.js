// API logic

import { Router } from "itty-router";
import { supabase, hash, t } from "./utils";

// API Router
export const routerAPI = Router({ base: "/api/v1" });

// -----------------------------------------------------------------------------
// Routes
// -----------------------------------------------------------------------------

// Analytics on tutorial progress
routerAPI.post("/ping", async request => {
	const data = await request.json();

	await supabase.from(t("pings")).insert([{
		ip: await hash(request.headers.get("CF-Connecting-IP")),
		tutorial: data.tutorial,
		step_from: data.from,
		step_to: data.to,
		playground: data.playground,
		example: data.example
	}]);

	return new Response("pong");
});

// -----------------------------------------------------------------------------
// Otherwise, 404
// -----------------------------------------------------------------------------

routerAPI.all("*", () => new Response("Not Found.", { status: 404 }));
