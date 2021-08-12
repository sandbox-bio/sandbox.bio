import { get } from "svelte/store";
import { CLI } from "../src/terminal/cli";
const $CLI = get(CLI);
let observed;

describe("Test variables", () => {
	before(async () => {
		console.log("Initializing Aioli");
		await $CLI.init({
			tools: ["samtools/1.10"]
		});
	});

	it("Set/Read variable", async () => {
		observed = await $CLI.exec(`abc=123`);
		observed = await $CLI.exec(`echo $abc`);
		expect(observed).to.equal("123");

		observed = await $CLI.exec(`echo $doesntexist`);
		expect(observed).to.equal("");

		observed = await $CLI.exec(`def=456`);
		observed = await $CLI.exec(`env`);
		expect(observed).to.equal("abc=123\ndef=456");
	});

	it("Concatenate variables", async () => {
		observed = await $CLI.exec(`abc=123`);
		observed = await $CLI.exec(`def=456`);
		observed = await $CLI.exec(`echo "this  $abc is a $def test"`);
		expect(observed).to.equal("this  123 is a 456 test");
	});

	it("Redirect to filename stored in variable", async () => {
		observed = await $CLI.exec(`abc=test.txt`);
		observed = await $CLI.exec(`echo 789 > $abc`);
		observed = await $CLI.exec(`cat $abc`);
		expect(observed).to.equal("789");
	});
});
