import path from "path";
import { sveltekit } from "@sveltejs/kit/vite";
import { defineConfig } from "vitest/config";

export default defineConfig({
	plugins: [sveltekit()],
	resolve: {
		alias: {
			$components: path.resolve(__dirname, "./src/components"),
			$stores: path.resolve(__dirname, "./src/stores"),
			$content: path.resolve(__dirname, "./src/content"),
			$thirdparty: path.resolve(__dirname, "./src/thirdparty"),
		}
	},
	test: {
		include: ["src/**/*.{test,spec}.{js,ts}"]
	},
	ssr: {
		// Avoids PopperJS error: "cannot use import statement outside a module" error
		noExternal: ["@popperjs/core"]
	}
});
