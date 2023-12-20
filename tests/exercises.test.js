// =============================================================================
// Test exercise functionality
// =============================================================================

import { test, expect } from "@playwright/test";
import { goToTutorial } from "./utils";

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
	await expect(page.locator(".bi-check-circle-fill")).toHaveCount(2);
	await expect(page.locator(".bi-circle")).toHaveCount(0);
});
