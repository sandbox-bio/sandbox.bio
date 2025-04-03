// -----------------------------------------------------------------------------
// Utilities
// -----------------------------------------------------------------------------

// MD5 hash (source: https://stackoverflow.com/a/64795218)
export async function hash(message) {
	if (!message) return;

	const msgUint8 = new TextEncoder().encode(message); // encode as (utf-8) Uint8Array
	const hashBuffer = await crypto.subtle.digest("MD5", msgUint8); // hash the message
	const hashArray = Array.from(new Uint8Array(hashBuffer)); // convert buffer to byte array
	const hashHex = hashArray.map((b) => b.toString(16).padStart(2, "0")).join(""); // convert bytes to hex string

	return hashHex;
}
