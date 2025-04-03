import path from "path";
import { sveltekit } from "@sveltejs/kit/vite";
import { defineConfig } from "vitest/config";

export default defineConfig({
	plugins: [sveltekit()],
	resolve: {
		alias: {
			$src: path.resolve(__dirname, "./src/"),
			$components: path.resolve(__dirname, "./src/components"),
			$stores: path.resolve(__dirname, "./src/stores"),
			$content: path.resolve(__dirname, "./src/content"),
			$routes: path.resolve(__dirname, "./src/routes"),
			$thirdparty: path.resolve(__dirname, "./src/thirdparty")
		}
	},
	test: {
		include: ["src/**/*.{test,spec}.{js,ts}"]
	}
});
