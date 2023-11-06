import { test } from "@playwright/test";
import { expectXterm, goToTerminal } from "./utils";

test("Basic commands", async ({ page }) => {
	await goToTerminal(page);

	await expectXterm(page, "pwd", "/root/tutorial");
	await expectXterm(page, "hostname", "localhost");
});

test("Interactive commands (man, less, vim, nano)", async ({ page }) => {
	await goToTerminal(page);

	await expectXterm(page, "man grep", "GREP(1)", ({ keyboard }) => keyboard.type("q"));
	await expectXterm(page, "less /root/.bashrc", "executed by bash(1) for non-login shells", ({ keyboard }) => keyboard.type("q"));
	await expectXterm(page, "nano /root/.bashrc", "GNU nano 7.2", ({ keyboard }) => keyboard.press("Control+X"));
	await expectXterm(page, "vim", "VIM - Vi IMproved", ({ keyboard }) => {
		keyboard.press(":")
		keyboard.type("q")
		keyboard.press("Enter")
	});
});
