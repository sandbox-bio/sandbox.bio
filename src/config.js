import { writable } from "svelte/store";
import { config as bedtoolsIntro } from "tutorials/bedtools-intro/config.js";
import { config as bowtie2Intro } from "tutorials/bowtie2-intro/config.js";

export const config = writable({
	status: "ready",
	tutorials: [
		bedtoolsIntro,
		bowtie2Intro
	]
});
