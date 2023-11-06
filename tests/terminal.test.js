// =============================================================================
// Test terminal commands and make sure tools installed in the v86 Dockerfile
// are accessible. These tests only check that the tool runs, not that the output
// is valid. e.g. a tool's --version might run but running it with another command
// might fail due to a compilation issue.
// =============================================================================

import { test } from "@playwright/test";
import { expectXterm, goToTerminal } from "./utils";

let page;
const tools = [
	// sandbox.bio v1: Installed with apt
	{ name: "jq", expected: "jq-1.6" },
	{ name: "tree", expected: "tree v2.1.0" },
	{ name: "seqtk", command: "seqtk", expected: "Version: 1.3-r106" },
	{ name: "fasttree", command: "fasttree 2>&1 | head", expected: "FastTree version 2.1.11" },
	{ name: "kalign", expected: "kalign 3.3.5" },
	{ name: "bedtools", expected: "bedtools v2.30.0" },
	{ name: "samtools", expected: "samtools 1.16.1" },
	{ name: "bcftools", expected: "bcftools 1.16" },
	{ name: "ivar", command: "ivar version", expected: "iVar version 1.3.1" },

	// sandbox.bio v1: Built from source
	{ name: "htsfile", expected: "htsfile (htslib) 1.18" },
	{ name: "bgzip", expected: "bgzip (htslib) 1.18" },
	{ name: "tabix", expected: "tabix (htslib) 1.18" },
	{ name: "fastp", expected: "fastp 0.20.1" },
	// { name: "lsd2", command: "lsd2 -V", expected: "lsd2 v.2.4.1" },
	// { name: "tn93", expected: "v1.0.11" },
	// { name: "bowtie2", expected: "bowtie2-align-s version 2.5.1" },
	{ name: "minimap2", expected: "2.26-r1175" },

	// sandbox.bio v2: Installed with apt
	{ name: "jellyfish", expected: "jellyfish 2.3.0" },
	{ name: "seqkit", command: "seqkit version", expected: "seqkit v2.3.0" },
	{ name: "freebayes", expected: "version:  v1.3.6" },
	{ name: "kraken2", expected: "Kraken version 2.1.2" },
	{ name: "lumpy", command: "lumpy 2>&1 | head", expected: "v 0.2.13" },
	{ name: "nucmer", expected: "NUCmer (NUCleotide MUMmer) version 3.1" },
	{ name: "sniffles", expected: "Sniffles2, Version 2.0.7" },
	{ name: "bioawk", expected: "awk version 20110810" },
	{ name: "vcflib", expected: "vcflength" },

	// sandbox.bio v2: Built from source
	{ name: "ViralMSA.py", command: "ViralMSA.py -h", expected: "current version: 1.1.32" },
	{ name: "csvtk", command: "csvtk | head", expected: "Version: 0.27.2" },
	{ name: "kallisto", expected: "kallisto 0.48.0" },
	{ name: "mmseqs", command: "mmseqs | head", expected: "MMseqs2 Version: 7e2840992948ee89dcc336522dc98a74fe0adf00" },
	{ name: "foldseek", command: "foldseek | head", expected: "foldseek Version: 946841ff3b15531349a9883358b3a3052b368da9" },
	{ name: "viral_consensus", expected: "viral_consensus v0.0.4" }
];

// Initialize terminal first
test.beforeAll(async ({ browser }) => {
	page = await browser.newPage();
	await goToTerminal(page);
  });
  
// Check that basic terminal commands work
test("Basic commands", async () => {
	await goToTerminal(page);

	await expectXterm(page, "pwd", "/root/tutorial");
	await expectXterm(page, "hostname", "localhost");
});

// Check that interactive commands work
test("Interactive commands (man, less, vim, nano)", async () => {
	await goToTerminal(page);

	await expectXterm(page, "man grep", "GREP(1)", ({ keyboard }) => keyboard.type("q"));
	await expectXterm(page, "less /root/.bashrc", "executed by bash(1) for non-login shells", ({ keyboard }) => keyboard.type("q"));
	await expectXterm(page, "nano /root/.bashrc", "GNU nano 7.2", ({ keyboard }) => keyboard.press("Control+X"));
	await expectXterm(page, "vim", "VIM - Vi IMproved", async ({ keyboard }) => {
		await keyboard.press(":");
		await keyboard.type("q");
		await keyboard.press("Enter");
	});
});

// Check that installed tools work
for (const tool of tools) {
	const command = tool.command || `${tool.name} --version`;
	test(`sandbox.bio command: ${command}`, async () => {
		await expectXterm(page, command, tool.expected);
	});
}
