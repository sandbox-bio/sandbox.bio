// =============================================================================
// Test exercise and quiz functionality
// =============================================================================

import { test, expect } from "@playwright/test";
import { goToTutorial } from "./utils";

test("Quiz validation - 1/2", async ({ browser }) => {
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
});

test("Quiz validation - 2/2", async ({ browser }) => {
	const page = await browser.newPage();
	await goToTutorial(page, "ifb-linux-basics-1", 2);

	// Make sure nothing is select in quiz initially
	await expect(page.locator(".alert-primary")).toHaveCount(1);

	// Choose wrong answer
	await page.getByLabel("No").check();
	await page.getByRole("button", { name: "Submit" }).click();

	// Should inform user
	await expect(page.locator(".alert-primary")).toHaveCount(1);
	await expect(page.locator(".alert-success")).toHaveCount(0);
	await page.waitForSelector(`:has-text("That doesn't look right")`);
});

test("Exercise validation", async ({ browser }) => {
	const page = await browser.newPage();
	await goToTutorial(page, "bedtools-intro", 14);

	// Initially, there should be no checkmarks
	await expect(page.locator(".bi-circle")).toHaveCount(2);
	await expect(page.locator(".bi-check-circle-fill")).toHaveCount(0);

	// Solve 1/2 of the criteria
	await page.keyboard.type("touch notexons.bed");
	await page.keyboard.press("Enter");
	await expect(page.locator(".bi-circle")).toHaveCount(1);
	await expect(page.locator(".bi-check-circle-fill")).toHaveCount(1);

	// Undo
	await page.keyboard.type("rm notexons.bed");
	await page.keyboard.press("Enter");
	await expect(page.locator(".bi-circle")).toHaveCount(2);
	await expect(page.locator(".bi-check-circle-fill")).toHaveCount(0);

	// Solve 2/2 of the criteria
	await page.keyboard.type("bedtools complement -i exons.bed -g genome.txt > notexons.bed");
	await page.keyboard.press("Enter");
	await expect(page.locator(".bi-circle")).toHaveCount(0);
	await expect(page.locator(".bi-check-circle-fill")).toHaveCount(2);
});
