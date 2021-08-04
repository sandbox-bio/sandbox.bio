import { readable } from "svelte/store";
import { config as bedtoolsIntro } from "tutorials/bedtools-intro/config.js";
import { config as bowtie2Intro } from "tutorials/bowtie2-intro/config.js";
import { config as samtoolsIntro } from "tutorials/samtools-intro/config.js";
import { config as hiddenMessage } from "tutorials/hidden-message/config.js";

export const config = readable({
	tutorials: [
		bedtoolsIntro,
		bowtie2Intro,
		samtoolsIntro,
		// hiddenMessage  // uncomment this + add "bcftools/1.10" to src/Terminal.svelte
	]
});
