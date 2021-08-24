// Test variable support

import { get } from "svelte/store";
import { CLI } from "../src/terminal/cli";
const $CLI = get(CLI);
let observed;
let expected;

describe("Test globbing", () => {
	before(async () => {
		console.log("Initializing Aioli");
		await $CLI.init({
			tools: ["samtools/1.10"]
		});
		await $CLI.exec(`echo abc > test1.txt; echo abc > test2.txt; echo abc > test3.bed`);
		await $CLI.exec(`mkdir folder1 folder2 folder3`);
		await $CLI.exec(`echo abc > folder1/file1-inside-folder1; echo abc > folder1/file2-inside-folder1`);
		await $CLI.exec(`echo abc > folder2/file1-inside-folder2; echo abc > folder2/file2-inside-folder2`);
		await $CLI.exec(`echo abc > folder3/file1-inside-folder3; echo abc > folder3/file2-inside-folder3`);
	});

	it("ls *", async () => {
		observed = await $CLI.exec(`ls *`);
		expected = await $CLI.exec(`ls`);
		expect(observed).to.equal(expected);
	});

	it("ls folder/*", async () => {
		observed = await $CLI.exec(`ls folder1/*`);
		expected = await $CLI.exec(`ls folder1/`);
		expect(observed).to.equal(expected);
	});

	it("ls f*", async () => {
		observed = await $CLI.exec(`ls f*`);
		expected = await $CLI.exec(`ls folder1 folder2 folder3`);
		expect(observed).to.equal(expected);
	});

	// Folders with patterns at the end can be tricky because of the implied "/"
	it("ls fold*3", async () => {
		observed = await $CLI.exec(`ls fol*3`);
		expected = await $CLI.exec(`ls folder3`);
		expect(observed).to.equal(expected);

		observed = await $CLI.exec(`ls fol*3/`);
		expected = await $CLI.exec(`ls folder3/`);
		expect(observed).to.equal(expected);
	});

	it("ls *pattern", async () => {
		observed = await $CLI.exec(`ls *.txt`);
		expect(observed).to.include("test1.txt");
		expect(observed).to.include("test2.txt");
		expect(observed).to.not.include("test3.bed");

		observed = await $CLI.exec(`ls *.bed`);
		expect(observed).to.not.include("test1.txt");
		expect(observed).to.not.include("test2.txt");
		expect(observed).to.include("test3.bed");

		observed = await $CLI.exec(`ls test*`);
		expect(observed).to.include("test1.txt");
		expect(observed).to.include("test2.txt");
		expect(observed).to.include("test3.bed");
	});

	it("ls ?pattern", async () => {
		observed = await $CLI.exec(`ls test?.txt`);
		expect(observed).to.include("test1.txt");
		expect(observed).to.include("test2.txt");
		expect(observed).to.not.include("test3.bed");
	});

	it("ls ?*pattern", async () => {
		observed = await $CLI.exec(`ls test?.*`);
		expect(observed).to.include("test1.txt");
		expect(observed).to.include("test2.txt");
		expect(observed).to.include("test3.bed");

		observed = await $CLI.exec(`ls test*.t?t`);
		expect(observed).to.include("test1.txt");
		expect(observed).to.include("test2.txt");
		expect(observed).to.not.include("test3.bed");
	});

	it("ls folder/*pattern", async () => {
		observed = await $CLI.exec(`ls folder1/*inside*`);
		expect(observed).to.include("file1-inside-folder1");
		expect(observed).to.include("file2-inside-folder1");
		expect(observed).to.not.include("file1-inside-folder2");
		expect(observed).to.not.include("file2-inside-folder2");
	});

	it("ls doesn*exist", async () => {
		// This should fail because the path doesn't exist
		try {
			await $CLI.exec(`ls test4*`)			
		} catch (error) {
			return;
		}
		expect(1).to.equal(2);
	});
});
