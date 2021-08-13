import { get } from "svelte/store";
import { CLI } from "../src/terminal/cli";

const $CLI = get(CLI);
const FILE_SAM = "/samtools/examples/toy.sam";

let observed;

describe("Test piping", () => {
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

});
