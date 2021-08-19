// API logic

import { Router } from "itty-router";

// API Router
export const routerAPI = Router({ base: "/api/v1" });

// Routes
routerAPI.get('/todos', () => new Response("Todos"));
routerAPI.get('/todos/:id', ({ params }) => new Response(`Todo #${params.id}`));
routerAPI.all('*', () => new Response('Not Found.', { status: 404 }));
