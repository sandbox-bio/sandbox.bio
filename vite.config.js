import { sveltekit } from "@sveltejs/kit/vite";
import { defineConfig } from "vitest/config";

export default defineConfig({
	plugins: [sveltekit()],
	test: {
		include: ["src/**/*.{test,spec}.{js,ts}"]
	},
	ssr: {
		// Avoids PopperJS error: "cannot use import statement outside a module" error
		noExternal: ["@popperjs/core"]
	}
});
