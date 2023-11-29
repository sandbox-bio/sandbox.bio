export async function handleError({ error, event }) {
	console.error("--------------------");
	console.error(error);
	console.error("--------------------");

	return {};
}
