// API logic

import { Router } from "itty-router";

// API Router
export const routerAPI = Router({ base: "/api/v1" });

// Routes
routerAPI.get('/todos', () => new Response('Todoss Index!'));
routerAPI.get('/todos/:id', ({ params }) => new Response(`Todo #${params.id}`));
routerAPI.post('/todos', async request => {
	const content = await request.json();
	return new Response('Creating Todo: ' + JSON.stringify(content));
})
routerAPI.all('*', () => new Response('Not Found.', { status: 404 }));
