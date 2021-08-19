module.exports = {
	target: "webworker",
	entry: "./index.js",
	mode: "production",
	// supabase.js relies on "cross-fetch", which uses XMLHttpRequest, which isn't supported
	externals: [{
		"cross-fetch": "fetch"
	}]
}
