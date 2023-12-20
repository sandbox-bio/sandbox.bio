// =============================================================================
// Test quiz functionality
// =============================================================================

import { test, expect } from "@playwright/test";
import { goToTutorial } from "./utils";

test("Quiz validation - One valid answer - Select correct answer", async ({ browser }) => {
	const page = await browser.newPage();

	await goToTutorial(page, "ifb-linux-basics-1", 2);

	// Make sure nothing is select in quiz initially
	await expect(page.locator(".alert-primary")).toHaveCount(1);

	// Choose right answer
	await page.getByLabel("Yes").check();
	await page.getByRole("button", { name: "Submit" }).click();

	// Should inform the user the answer is correct
	await expect(page.locator(".alert-primary")).toHaveCount(0);
	await expect(page.locator(".alert-success")).toHaveCount(1);
	await page.waitForSelector(`:has-text("That is correct!")`);

	// Refresh the page to make sure we save state locally
	await goToTutorial(page, "ifb-linux-basics-1", 2);
	await expect(page.locator(".alert-primary")).toHaveCount(0);
	await expect(page.locator(".alert-success")).toHaveCount(1);
	await page.waitForSelector(`:has-text("That is correct!")`);
});

test("Quiz validation - One valid answer - Select wrong answer", async ({ browser }) => {
	const page = await browser.newPage();

	// User provides the wrong answer
	await goToTutorial(page, "ifb-linux-basics-1", 3);

	// Make sure nothing is select in quiz initially
	await expect(page.locator(".alert-primary")).toHaveCount(1);

	// Choose wrong answer
	await page.getByLabel("Terminal", { exact: true }).check();
	await page.getByRole("button", { name: "Submit" }).click();

	// Should inform user
	await expect(page.locator(".alert-primary")).toHaveCount(1);
	await expect(page.locator(".alert-success")).toHaveCount(0);
	await page.waitForSelector(`:has-text("That doesn't look right")`);
});

test("Quiz validation - Multiple valid answers - Select correct answer", async ({ browser }) => {
	const page = await browser.newPage();

	// User provides the correct answer
	await goToTutorial(page, "ifb-linux-basics-2", 2);

	// Make sure nothing is select in quiz initially
	await expect(page.locator(".alert-primary")).toHaveCount(1);

	// Choose right answer
	await page.getByLabel("hg19").check();
	await page.getByLabel("hg38").check();
	await page.getByRole("button", { name: "Submit" }).click();

	// Should inform the user the answer is correct
	await expect(page.locator(".alert-primary")).toHaveCount(0);
	await expect(page.locator(".alert-success")).toHaveCount(1);
	await page.waitForSelector(`:has-text("That is correct!")`);

	// Refresh the page to make sure we save state locally
	await goToTutorial(page, "ifb-linux-basics-2", 2);
	await expect(page.locator(".alert-primary")).toHaveCount(0);
	await expect(page.locator(".alert-success")).toHaveCount(1);
	await page.waitForSelector(`:has-text("That is correct!")`);
});

test("Quiz validation - Multiple valid answers - Select wrong answer", async ({ browser }) => {
	const page = await browser.newPage();

	await goToTutorial(page, "ifb-linux-basics-3", 7);

	// Make sure nothing is select in quiz initially
	await expect(page.locator(".alert-primary")).toHaveCount(1);

	// Choose wrong answer
	await page.getByLabel("cut -f 2,3,6 MACS2.csv").check();
	await page.getByRole("button", { name: "Submit" }).click();

	// Should inform user
	await expect(page.locator(".alert-primary")).toHaveCount(1);
	await expect(page.locator(".alert-success")).toHaveCount(0);
	await page.waitForSelector(`:has-text("That doesn't look right")`);
});
