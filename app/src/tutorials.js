import { readable } from "svelte/store";
import { config as bedtoolsIntro } from "tutorials/bedtools-intro/config.js";
import { config as bowtie2Intro } from "tutorials/bowtie2-intro/config.js";
import { config as samtoolsIntro } from "tutorials/samtools-intro/config.js";
import { config as dnaSecrets } from "tutorials/dna-secrets/config.js";

export const tutorials = readable([
	bedtoolsIntro,
	bowtie2Intro,
	samtoolsIntro,
	dnaSecrets
]);
