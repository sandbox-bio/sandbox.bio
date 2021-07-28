import Aioli from "@biowasm/aioli";
import { CoreUtils } from "../../src/terminal/coreutils"

const FILE_SAM = "/samtools/examples/toy.sam";
const FILE_FA = "/samtools/examples/toy.fa";
let observed, expected;

describe("Test coreutils", () => {
	before(async () => {
		console.log("Initializing Aioli");
		CoreUtils.CLI = await new Aioli("samtools/1.10", { env: "stg", debug: true });
	});

	it("ls (folder)", async () => {
		observed = (await CoreUtils.ls(["/"], true)).map(d => d.name);
		expected = ["tmp/", "home/", "dev/", "proc/", "samtools/", "shared/"];
		expect(JSON.stringify(observed)).to.equal(JSON.stringify(expected));

		observed = (await CoreUtils.ls([FILE_SAM], true)).map(d => d.name);
		expected = ["toy.sam"];
		expect(JSON.stringify(observed)).to.equal(JSON.stringify(expected));
	});

	it("head", async () => {
		observed = (await CoreUtils.head([FILE_SAM])).split("\n").length;
		expected = 10;
		expect(observed).to.equal(expected);

		observed = (await CoreUtils.head(["-n", 5, FILE_SAM])).split("\n").length;
		expected = 5;
		expect(observed).to.equal(expected);
	});

	it("tail", async () => {
		observed = (await CoreUtils.tail([FILE_SAM])).split("\n").length;
		expected = 11;
		expect(observed).to.equal(expected);

		observed = (await CoreUtils.tail(["-n", 5, FILE_SAM])).split("\n").length;
		expected = 6;
		expect(observed).to.equal(expected);
	});

	it("wc", async () => {
		observed = await CoreUtils.wc([FILE_SAM]);
		expected = 786;
		expect(observed).to.equal(expected);

		observed = await CoreUtils.wc(["-c", FILE_SAM]);
		expected = 786;
		expect(observed).to.equal(expected);

		observed = await CoreUtils.wc(["-l", FILE_SAM]);
		expected = 15;
		expect(observed).to.equal(expected);
	});

	it("cat", async () => {
		observed = await CoreUtils.cat([FILE_FA]);
		expected = `>ref\nAGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT\n>ref2\naggttttataaaacaattaagtctacagagcaactacgcg\n`;
		expect(observed).to.equal(expected);
	});

	it("cd / pwd", async () => {
		observed = await CoreUtils.pwd();
		expected = "/shared/data";
		expect(observed).to.equal(expected);

		// Change directory
		await CoreUtils.cd([ "/tmp" ]);
		observed = await CoreUtils.pwd();
		expected = "/tmp";
		expect(observed).to.equal(expected);
		await CoreUtils.cd([ "/shared/data" ]);
	});

	it("echo", async () => {
		observed = await CoreUtils.echo(["output something"]);
		expected = `output something`;
		expect(observed).to.equal(expected);
	});

	it("mkdir / rmdir", async () => {
		await CoreUtils.mkdir(["a", "b", "c"]);
		observed = (await CoreUtils.ls(["."], true)).map(d => d.name);
		expected = ["a/", "b/", "c/"];
		expect(JSON.stringify(observed)).to.equal(JSON.stringify(expected));

		await CoreUtils.rmdir(["a", "b"]);
		observed = (await CoreUtils.ls(["."], true)).map(d => d.name);
		expected = ["c/"];
		expect(JSON.stringify(observed)).to.equal(JSON.stringify(expected));
	});

	it("mv / rm", async () => {
		observed = (await CoreUtils.ls([FILE_SAM], true)).map(d => d.name);
		expected = ["toy.sam"];
		expect(JSON.stringify(observed)).to.equal(JSON.stringify(expected));

		// There should be a new file once we move it
		await CoreUtils.mv([FILE_SAM, FILE_SAM + ".2"]);
		observed = (await CoreUtils.ls([FILE_SAM + ".2"], true)).map(d => d.name);
		expected = ["toy.sam.2"];
		expect(JSON.stringify(observed)).to.equal(JSON.stringify(expected));

		// And the old one shouldn't be there anymore
		try {
			observed = (await CoreUtils.ls([FILE_SAM], true)).map(d => d.name);
		} catch (error) {
			expect(error).to.equal(`${FILE_SAM}: No such file or directory`);
		}

		// Delete the file
		await CoreUtils.rm([FILE_SAM + ".2"]);
		try {
			observed = (await CoreUtils.ls([FILE_SAM + ".2"], true)).map(d => d.name);
		} catch (error) {
			expect(error).to.equal(`${FILE_SAM + ".2"}: No such file or directory`);
		}
	});
});
