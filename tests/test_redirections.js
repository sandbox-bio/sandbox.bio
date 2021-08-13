import { get } from "svelte/store";
import { CLI } from "../src/terminal/cli";

const $CLI = get(CLI);
const FILE_SAM = "/samtools/examples/toy.sam";

let observed;

describe("Test redirections", () => {
	before(async () => {
		console.log("Initializing Aioli");
		await $CLI.init({
			tools: ["samtools/1.10"]
		});
	});

	it("Single pipe", async () => {
		observed = await $CLI.exec(`samtools view | head -n 1`);
		expect(observed).to.equal("Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]");
	});

	it("Single pipe + file redirection", async () => {
		observed = await $CLI.exec(`samtools view | head -n 1 > tmp`);
		observed = await $CLI.exec(`cat tmp`);
		expect(observed).to.equal("Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]");
	});

	it("Multiple pipes", async () => {
		observed = await $CLI.exec(`samtools view | head -n 7 | wc -l`);
		expect(observed).to.equal("7");
	});

	it("Multiple pipes + file redirection + output is not a string", async () => {
		observed = await $CLI.exec(`samtools view | head -n 7 | wc -l > tmp`);
		observed = await $CLI.exec(`cat tmp`);
		expect(observed).to.equal("7");
	});

	it("Process substitution", async () => {
		observed = await $CLI.exec(`wc -l <(head -n 7 ${FILE_SAM})`);
		expect(observed).to.equal("7");
	});

	it("Command substitution", async () => {
		observed = await $CLI.exec(`head -n $(head -n 7 ${FILE_SAM} | wc -l) ${FILE_SAM} | wc -l`);
		expect(observed).to.equal("7");
	});

	it("Write/Append to file", async () => {
		observed = await $CLI.exec(`echo "test" > tmp`);
		observed = await $CLI.exec(`echo "test" > tmp`);
		observed = await $CLI.exec(`cat tmp`);
		expect(observed).to.equal("test");

		observed = await $CLI.exec(`echo "test" > tmp`);
		observed = await $CLI.exec(`echo "test" >> tmp`);
		observed = await $CLI.exec(`cat tmp`);
		expect(observed).to.equal("test\ntest");
	});

	it("Using '-' as a replacement for stdin", async () => {
		observed = await $CLI.exec(`echo "test" | cat -`);
		expect(observed).to.equal("test");

		observed = await $CLI.exec(`echo "test" | cat`);
		expect(observed).to.equal("test");

		observed = await $CLI.exec(`samtools view /samtools/examples/toy.sam | cat - | wc -l`);
		expect(observed).to.equal("12");
	});

	it("Treating stderr correctly", async () => {
		let stderr = "";
		const samtoolsHelp = `\nProgram: samtools (Tools for alignments in the SAM format)\nVersion: 1.10 (using htslib 1.10)\n\nUsage:   samtools <command> [options]\n\nCommands:\n  -- Indexing\n     dict           create a sequence dictionary file\n     faidx          index/extract FASTA\n     fqidx          index/extract FASTQ\n     index          index alignment\n\n  -- Editing\n     calmd          recalculate MD/NM tags and \'=\' bases\n     fixmate        fix mate information\n     reheader       replace BAM header\n     targetcut      cut fosmid regions (for fosmid pool only)\n     addreplacerg   adds or replaces RG tags\n     markdup        mark duplicates\n\n  -- File operations\n     collate        shuffle and group alignments by name\n     cat            concatenate BAMs\n     merge          merge sorted alignments\n     mpileup        multi-way pileup\n     sort           sort alignment file\n     split          splits a file by read group\n     quickcheck     quickly check if SAM/BAM/CRAM file appears intact\n     fastq          converts a BAM to a FASTQ\n     fasta          converts a BAM to a FASTA\n\n  -- Statistics\n     bedcov         read depth per BED region\n     coverage       alignment depth and percent coverage\n     depth          compute the depth\n     flagstat       simple stats\n     idxstats       BAM index stats\n     phase          phase heterozygotes\n     stats          generate stats (former bamcheck)\n\n  -- Viewing\n     flags          explain BAM flags\n     tview          text alignment viewer\n     view           SAM<->BAM<->CRAM conversion\n     depad          convert padded BAM to unpadded BAM\n\n`;

		// samtools outputs its main help screen to stderr
		observed = await $CLI.exec("samtools", d => stderr = d);
		expect(observed).to.equal("");
		expect(stderr).to.equal(samtoolsHelp);

		observed = await $CLI.exec("samtools | head -n 1", d => stderr = d);
		expect(observed).to.equal("");
		expect(stderr).to.equal(samtoolsHelp);

		observed = await $CLI.exec("samtools > somefile", d => stderr = d);
		expect(observed).to.equal("");
		expect(stderr).to.equal(samtoolsHelp);

		// samtools --version-only outputs to stdout
		observed = await $CLI.exec("samtools --version-only");
		expect(observed).to.equal("1.10+htslib-1.10");
	});

});
