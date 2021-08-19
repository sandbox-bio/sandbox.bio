// API logic

import { Router } from "itty-router";
import { createClient } from "@supabase/supabase-js";

// API Router
export const routerAPI = Router({ base: "/api/v1" });

// Supabase Client
const supabase = createClient(SUPABASE_URL, SUPABASE_API_KEY);  // stored in Cloudflare KV secrets

// Routes
routerAPI.get('/todos', () => new Response("Todos"));
routerAPI.get('/todos/:id', ({ params }) => new Response(`Todo #${params.id}`));
routerAPI.get('/tutorials', async request => {
	const { data, error } = await supabase.from("tutorials").select();
	console.log("sup")

	return new Response(JSON.stringify(data));
})
routerAPI.all('*', () => new Response('Not Found.', { status: 404 }));
