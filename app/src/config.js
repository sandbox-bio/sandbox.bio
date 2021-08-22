import { writable } from "svelte/store";
import { config as bedtoolsIntro } from "tutorials/bedtools-intro/config.js";
import { config as bowtie2Intro } from "tutorials/bowtie2-intro/config.js";
import { config as samtoolsIntro } from "tutorials/samtools-intro/config.js";
import { config as dnaSecrets } from "tutorials/dna-secrets/config.js";

const URL_API = `https://${window.location.hostname != "localhost" ? window.location.hostname : "dev.sandbox.bio"}/api/v1`;

export const config = writable({
	tutorials: [
		bedtoolsIntro,
		bowtie2Intro,
		samtoolsIntro,
		dnaSecrets
	],
	api: URL_API,
	cli: {
		user: "guest",
		hostname: "sandbox"
	}
});
