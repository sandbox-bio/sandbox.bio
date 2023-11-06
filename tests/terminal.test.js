import { test } from "@playwright/test";
import { expectXterm, goToTerminal } from "./utils";

test("Basic terminal commands", async ({ page }) => {
	await goToTerminal(page);

	await expectXterm(page, "pwd", "/root/tutorial");
	await expectXterm(page, "man grep", "GREP(1)", ({ xterm }) => xterm.press("q"));

	// Make sure the expectXterm callback above worked
	await expectXterm(page, "hostname", "localhost");
});
