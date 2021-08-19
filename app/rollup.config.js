import svelte from "rollup-plugin-svelte";
import commonjs from "@rollup/plugin-commonjs";
import resolve from "@rollup/plugin-node-resolve";
import livereload from "rollup-plugin-livereload";
import includePaths from "rollup-plugin-includepaths";
import css from "rollup-plugin-css-only";
import { terser } from "rollup-plugin-terser";
import { markdown } from "svelte-preprocess-markdown";

const production = !process.env.ROLLUP_WATCH;

function serve() {
	let server;

	function toExit() {
		if (server) server.kill(0);
	}

	return {
		writeBundle() {
			if (server) return;
			server = require("child_process").spawn("npm", ["run", "start", "--", "--dev"], {
				stdio: ["ignore", "inherit", "inherit"],
				shell: true
			});

			process.on("SIGTERM", toExit);
			process.on("exit", toExit);
		}
	};
}

export default {
	input: "src/main.js",
	output: {
		sourcemap: true,
		format: "iife",
		name: "app",
		file: "public/build/bundle.js"
	},
	plugins: [
		svelte({
			// Note that exercises can be written in a mix of Svelte/Markdown
			extensions: [".svelte", ".md"],
			preprocess: markdown(),
			// Run-time checks when not in production
			compilerOptions: { dev: !production }
		}),

		// Extract component CSS into a separate file - better for performance
		css({ output: "bundle.css" }),

		// Resolve external dependencies installed from npm
		resolve({ browser: true, dedupe: ["svelte"] }),
		commonjs(),

		// Launch server in dev mode
		!production && serve(),
		!production && livereload("public"),

		// Minify
		production && terser(),

		// Define other include paths so the code is more terse
		includePaths({ paths: ["./src/"] })
	],
	watch: {
		clearScreen: false
	}
};
