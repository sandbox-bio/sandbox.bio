import { readable } from "svelte/store";
import { config as bedtoolsIntro } from "tutorials/bedtools-intro/config.js";
import { config as bowtie2Intro } from "tutorials/bowtie2-intro/config.js";
import { config as samtoolsIntro } from "tutorials/samtools-intro/config.js";

export const config = readable({
	tutorials: [
		bedtoolsIntro,
		bowtie2Intro,
		samtoolsIntro
	]
});
