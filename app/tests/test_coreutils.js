// Test coreutils

import { get } from "svelte/store";
import { CLI } from "../src/terminal/cli"; const $CLI = get(CLI);

const FILE_SAM = "/samtools/examples/toy.sam";
const FILE_FA = "/samtools/examples/toy.fa";
let observed;
let expected;

describe("Test coreutils", () => {
	before(async () => {
		console.log("Initializing Aioli");
		await $CLI.init({
			tools: ["samtools/1.10"]
		});
	});

	it("ls (folder)", async () => {
		observed = await $CLI.exec(`ls ${FILE_SAM}`);
		expect(observed).to.contain("toy.sam");
		expect(observed).to.contain("786 B");

		observed = await $CLI.exec("ls /");
		expect(observed).to.contain("tmp/");
		expect(observed).to.contain("home/");
	});

	it("head", async () => {
		observed = (await $CLI.exec(`head ${FILE_SAM}`)).split("\n").length;
		expect(observed).to.equal(10);

		observed = (await $CLI.exec(`head -n 5 ${FILE_SAM}`)).split("\n").length;
		expect(observed).to.equal(5);

		observed = await $CLI.exec(`head -n 1 ${FILE_SAM}`);
		expect(observed).to.equal("@SQ\tSN:ref\tLN:45");
	});

	it("tail", async () => {
		observed = (await $CLI.exec(`tail ${FILE_SAM}`)).split("\n").length;
		expect(observed).to.equal(10);

		observed = (await $CLI.exec(`tail -n 5 ${FILE_SAM}`)).split("\n").length;
		expect(observed).to.equal(5);

		observed = await $CLI.exec(`tail -n 1 ${FILE_SAM}`);
		expect(observed).to.equal("x6\t0\tref2\t14\t30\t23M\t*\t0\t0\tTaattaagtctacagagcaacta\t???????????????????????");
	});

	// FIXME: technically this should be 786 bytes, but modifying utils.readFiles messes up head/tail
	it("wc", async () => {
		observed = await $CLI.exec(`wc ${FILE_SAM}`);
		expect(observed).to.equal("14\t139\t785");

		observed = await $CLI.exec(`wc -l ${FILE_SAM}`);
		expect(observed).to.equal("14");

		observed = await $CLI.exec(`wc -w ${FILE_SAM}`);
		expect(observed).to.equal("139");

		observed = await $CLI.exec(`wc -c ${FILE_SAM}`);
		expect(observed).to.equal("785");
	});

	it("cat", async () => {
		observed = await $CLI.exec(`cat ${FILE_FA}`);
		expect(observed).to.equal(`>ref\nAGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT\n>ref2\naggttttataaaacaattaagtctacagagcaactacgcg`);
	});

	it("cd / pwd", async () => {
		observed = await $CLI.exec("pwd");
		expect(observed).to.equal("/shared/data");

		// Change directory
		await $CLI.exec("cd /tmp");
		observed = await $CLI.exec("pwd");
		expect(observed).to.equal("/tmp");

		// Test cd ~
		await $CLI.exec("cd ~");
		observed = await $CLI.exec("pwd");
		expect(observed).to.equal("/shared/data");

		// Test cd -
		await $CLI.exec("cd /shared");
		await $CLI.exec("cd /");
		await $CLI.exec("cd -");
		observed = await $CLI.exec("pwd");
		expect(observed).to.equal("/shared");
	});

	it("echo", async () => {
		observed = await $CLI.exec("echo 'output something'");
		expect(observed).to.equal("output something");
	});


	it("mktemp", async () => {
		observed = await $CLI.exec(`mktemp`);
		expect(observed).to.contain("/shared/tmp/tmp");

		observed = await $CLI.exec(`ls ${observed}`);  // throws error if file doesn't exist
		expect(observed).to.contain("tmp");
	});

	it("grep", async () => {
		observed = await $CLI.exec(`grep "Taa" ${FILE_SAM} | wc -l`);
		expect(observed).to.equal("4");

		observed = await $CLI.exec(`grep "TAA" ${FILE_SAM} | wc -l`);
		expect(observed).to.equal("3");

		observed = await $CLI.exec(`grep -i "TAA" ${FILE_SAM} | wc -l`);
		expect(observed).to.equal("9");

		observed = await $CLI.exec(`grep -i -v "TAA" ${FILE_SAM} | wc -l`);
		expect(observed).to.equal("5");
	});

	it("mkdir / rmdir", async () => {
		await $CLI.exec("mkdir a b c");

		observed = await $CLI.exec("ls");
		expect(observed).to.contain("a/");
		expect(observed).to.contain("b/");
		expect(observed).to.contain("c/");

		await $CLI.exec("rmdir a c");

		observed = await $CLI.exec("ls");
		expect(observed).to.contain("b/");
	});

	it("cp", async () => {
		// Copy file and make sure contents are identical
		await $CLI.exec(`cp ${FILE_SAM} copied.sam`);
		observed = await $CLI.exec(`cat copied.sam`);
		expected = await $CLI.exec(`cat ${FILE_SAM}`);
		expect(observed).to.equal(expected);

		// Overwrite and make sure it still works
		await $CLI.exec(`echo "hello" > somefile`);
		await $CLI.exec(`cp somefile ${FILE_SAM}`);
		observed = await $CLI.exec(`cat ${FILE_SAM}`);
		expect(observed).to.equal("hello");
	});

	// Run this test last since we're deleting files :)
	it("mv / rm", async () => {
		observed = await $CLI.exec(`ls ${FILE_SAM}`);
		expect(observed).to.contain("toy.sam");

		// There should be a new file once we move it
		const fileBlam = FILE_SAM.replace("toy.sam", "toy.blam");
		await $CLI.exec(`mv ${FILE_SAM} ${fileBlam}`);
		observed = await $CLI.exec(`ls /samtools/examples/`);
		expect(observed).to.contain("toy.blam");
		expect(observed).to.not.contain("toy.sam");

		// And the old one shouldn't be there anymore
		try {
			observed = await $CLI.exec(`ls ${FILE_SAM}`);
		} catch (error) {
			expect(error).to.equal(`${FILE_SAM}: No such file or directory`);
		}

		// Delete the file
		await $CLI.exec(`rm ${fileBlam}`)
		try {
			observed = await $CLI.exec(`ls ${fileBlam}`);
		} catch (error) {
			expect(error).to.equal(`${fileBlam}: No such file or directory`);
		}
	});
});
