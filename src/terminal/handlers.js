// =============================================================================
// xterm.js handlers
// =============================================================================

// Keyboard shortcuts
export function handleShortcuts(key)
{
	// Ctrl + L = Clear terminal
	if(key.domEvent.ctrlKey && key.domEvent.key == "l")
		term.write(`${ANSI_CLEAR}$ `);

	// Ctrl + A = Beginning of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "a")
		termEcho.setCursor(0);

	// Ctrl + E = End of line
	if(key.domEvent.ctrlKey && key.domEvent.key == "e")
		termEcho.setCursor(Infinity);
}

// Auto-completes common commands
export function handleAutocomplete(index, tokens)
{
	const command = tokens[0];
	console.log(index, tokens);

	// Root autocomplete
	if(index == 0)
		return ["samtools", "bedtools2"];

	// Samtools autocomplete. Need `&& tokens[index]`, otherwise results in "samtools samtools"
	if(index == 1 && command == "samtools" && tokens[index])
		return ["view", "index", "sort"];

	return [];
}
