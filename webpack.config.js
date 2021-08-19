module.exports = {
	target: "webworker",
	entry: "./index.js",
	mode: "production",
	externals: [{
		// supabase.js uses "cross-fetch", which uses "XMLHttpRequest", which
		// isn't supported in Cloudflare Workers.
		"cross-fetch": "fetch"
	}]
}
