// Initialize xterm.js library and associated addons

import { readable } from "svelte/store";
import { Terminal } from "xterm";
import { WebLinksAddon } from "xterm-addon-web-links";
import { SerializeAddon } from "xterm-addon-serialize";
import { FitAddon } from "xterm-addon-fit";
import LocalEchoController from "local-echo";


// =============================================================================
// Initialize xterm
// =============================================================================

// Xterm.js
const term = new Terminal({
	convertEol: true,
	cursorBlink: true
});

// Xterm.js add-ons
const addons = {
	echo: new LocalEchoController(null, {                  // Echo controller
		historySize: 1000
	}),
	serialize: new SerializeAddon(),                       // Can be used to save state
	links: new WebLinksAddon(),                            // Supports links in the terminal
	fit: new FitAddon()                                    // Makes terminal fit HTML element
};

// Attach addons
for(let addonName in addons)
	term.loadAddon(addons[addonName]);

// Disable local-echo's handling of tabs for autocompletion. See handleAutocomplete() for more info.
addons.echo.handleData_ = addons.echo.handleData;
addons.echo.handleData = (data) => {
	if(data == "\t")
		return;
	return addons.echo.handleData_(data);
}

// =============================================================================
// Export as readable stores
// =============================================================================

export const xterm = readable(term);
export const xtermAddons = readable(addons);
