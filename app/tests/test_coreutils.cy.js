// Test coreutils

import { get } from "svelte/store";
import { xtermAddons } from "../src/terminal/xterm";
import { CLI } from "../src/terminal/cli"; const $CLI = get(CLI);
import { TOOLS } from "./utils";

const FOLDER_SAMTOOLS = "/shared/samtools";
const FILE_SAM = `${FOLDER_SAMTOOLS}/examples/toy.sam`;
const FILE_FA = `${FOLDER_SAMTOOLS}/examples/toy.fa`;
const FILE_FA_CONTENTS = `>ref\nAGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT\n>ref2\naggttttataaaacaattaagtctacagagcaactacgcg\n`;
let observed;
let expected;

describe("Test coreutils", () => {
	before(async () => {
		console.log("Initializing Aioli");
		await $CLI.init({ tools: TOOLS });
	});

	it("ls (folder)", async () => {
		observed = await $CLI.exec(`ls -lah ${FILE_SAM}`);
		expect(observed).to.contain("toy.sam");
		expect(observed).to.contain("786");

		observed = await $CLI.exec("ls /");
		expect(observed).to.contain("tmp");
		expect(observed).to.contain("home");
	});

	it("head", async () => {
		observed = (await $CLI.exec(`head ${FILE_SAM}`)).trim().split("\n").length;
		expect(observed).to.equal(10);

		observed = (await $CLI.exec(`head -n 5 ${FILE_SAM}`)).trim().split("\n").length;
		expect(observed).to.equal(5);

		observed = await $CLI.exec(`head -n 1 ${FILE_SAM}`);
		expect(observed).to.equal("@SQ\tSN:ref\tLN:45\n");
	});

	it("tail", async () => {
		observed = (await $CLI.exec(`tail ${FILE_SAM}`)).trim().split("\n").length;
		expect(observed).to.equal(10);

		observed = (await $CLI.exec(`tail -n 5 ${FILE_SAM}`)).trim().split("\n").length;
		expect(observed).to.equal(5);

		observed = await $CLI.exec(`tail -n 1 ${FILE_SAM}`);
		expect(observed).to.equal("x6\t0\tref2\t14\t30\t23M\t*\t0\t0\tTaattaagtctacagagcaacta\t???????????????????????\n");
	});

	it("wc", async () => {
		observed = await $CLI.exec(`wc ${FILE_SAM}`);
		expect(observed).to.equal(` 14 139 786 ${FILE_SAM}\n`);

		observed = await $CLI.exec(`wc -l ${FILE_SAM}`);
		expect(observed).to.equal(`14 ${FILE_SAM}\n`);

		observed = await $CLI.exec(`wc -w ${FILE_SAM}`);
		expect(observed).to.equal(`139 ${FILE_SAM}\n`);

		observed = await $CLI.exec(`wc -c ${FILE_SAM}`);
		expect(observed).to.equal(`786 ${FILE_SAM}\n`);
	});

	it("cat", async () => {
		observed = await $CLI.exec(`cat ${FILE_FA}`);
		expect(observed).to.equal(FILE_FA_CONTENTS);
	});

	it("cd / pwd", async () => {
		observed = await $CLI.exec("pwd");
		expect(observed).to.equal("/shared/data\n");

		// Change directory
		await $CLI.exec("cd /tmp");
		observed = await $CLI.exec("pwd");
		expect(observed).to.equal("/tmp\n");

		// Test cd ~
		await $CLI.exec("HOME=/shared/data")
		await $CLI.exec("cd ~");
		observed = await $CLI.exec("pwd");
		expect(observed).to.equal("/shared/data\n");

		// Test ~ with other commands
		await $CLI.exec("cd /");
		await $CLI.exec(`cp ${FILE_FA} ~/test.fa`);
		observed = await $CLI.exec(`cat ~/test.fa`);
		expect(observed).to.equal(FILE_FA_CONTENTS);

		// Test cd -
		await $CLI.exec("cd /shared");
		await $CLI.exec("cd /");
		await $CLI.exec("cd -");
		observed = await $CLI.exec("pwd");
		expect(observed).to.equal("/shared\n");
	});

	it("echo", async () => {
		observed = await $CLI.exec("echo 'output something'");
		expect(observed).to.equal("output something\n");
	});

	it("mktemp", async () => {
		observed = await $CLI.exec(`mktemp`);
		expect(observed).to.contain("/shared/tmp/tmp");

		observed = await $CLI.exec(`ls ${observed}`);  // throws error if file doesn't exist
		expect(observed).to.contain("tmp");
	});

	it("grep+cut", async () => {
		observed = await $CLI.exec(`grep "Taa" ${FILE_SAM} | wc -l | cut -f1 -d' '`);
		expect(observed.trim()).to.equal("4");

		observed = await $CLI.exec(`grep "TAA" ${FILE_SAM} | wc -l | cut -f1 -d' '`);
		expect(observed.trim()).to.equal("3");

		observed = await $CLI.exec(`grep -i "TAA" ${FILE_SAM} | wc -l | cut -f1 -d' '`);
		expect(observed.trim()).to.equal("9");

		observed = await $CLI.exec(`grep -i -v "TAA" ${FILE_SAM} | wc -l | cut -f1 -d' '`);
		expect(observed.trim()).to.equal("5");
	});

	it("mkdir / rmdir", async () => {
		await $CLI.exec("mkdir a b c");
		observed = await $CLI.exec("ls");
		expect(observed).to.contain("a");
		expect(observed).to.contain("b");
		expect(observed).to.contain("c");

		await $CLI.exec("mkdir -p d/e/f g/h/i");
		observed = await $CLI.exec("ls -R d");
		expect(observed).to.contain("d:\ne\n\nd/e:\nf\n\nd/e/f:");
		observed = await $CLI.exec("ls -R g");
		expect(observed).to.contain("g:\nh\n\ng/h:\ni\n\ng/h/i:");

		await $CLI.exec("rmdir a c");
		observed = await $CLI.exec("ls");
		expect(observed).to.contain("b");
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
		expect(observed).to.equal("hello\n");
	});

	// Run this test last since we're deleting files :)
	it("mv / rm", async () => {
		observed = await $CLI.exec(`ls ${FILE_SAM}`);
		expect(observed).to.contain("toy.sam");

		// There should be a new file once we move it
		const fileBlam = FILE_SAM.replace("toy.sam", "toy.blam");
		await $CLI.exec(`mv ${FILE_SAM} ${fileBlam}`);
		observed = await $CLI.exec(`ls ${FOLDER_SAMTOOLS}/examples/`);
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

	it("whoami", async () => {
		await $CLI.exec("USER=someuser");
		observed = await $CLI.exec("whoami");
		expect(observed).to.equal("someuser\n");
	});

	it("touch", async () => {
		try {
			observed = await $CLI.exec("ls somefile");
			expect(1).to.equal(2);  // test fails if ls does not throw an error
		} catch (error) {}

		await $CLI.exec("touch somefile");
		observed = await $CLI.exec("ls somefile");  // test fails if ls throws an error
	});

	it("history", async () => {
		// Set history programmatically since it's controlled by the xterm.js UI
		const historyController = get(xtermAddons)?.echo?.history;
		historyController.entries = ["ls", "pwd", "hostname", "env", "history"];
		historyController.cursor = historyController.entries.length;

		observed = await $CLI.exec("history");
		expect(observed).to.equal("1\tls\n2\tpwd\n3\thostname\n4\tenv\n5\thistory\n");

		// Delete line 4 ("env")
		await $CLI.exec("history -d 4");
		observed = await $CLI.exec("history");
		expect(observed).to.equal("1\tls\n2\tpwd\n3\thostname\n4\thistory\n");

		// Clear history
		await $CLI.exec("history -c");
		observed = await $CLI.exec("history");
		expect(observed).to.equal("\n");
	});
});
