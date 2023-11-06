export async function goToTerminal(page) {
	// Open terminal and wait till it's ready
	await page.goto("/tutorials/playground");
	await page.waitForSelector(`:has-text("root@localhost")`);
}

export async function expectXterm(page, command, expect, callback) {
    // Type command one character at a time. Note that `fill` and
    // `pressSequentially` don't work for entering input in xterm.js
	await page.keyboard.type(command);
	await page.keyboard.press('Enter');

	// Validate
	await page.waitForSelector(`:has-text("${expect}")`);
	if (callback) {
		await callback({
			keyboard: page.keyboard,
		});
	}
}
