import adapter from "@sveltejs/adapter-vercel";
import { markdown } from "svelte-preprocess-markdown";

const config = {
	kit: {
		adapter: adapter()
	},
	extensions: [".svelte", ".md"],
	preprocess: [markdown()]
};

export default config;
