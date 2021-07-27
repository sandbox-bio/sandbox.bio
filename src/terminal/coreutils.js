// =============================================================================
// Coreutils
// =============================================================================

// Simulating a few basic utils
export class CoreUtils
{
	static CLI;

	// -------------------------------------------------------------------------
	// File system utils
	// -------------------------------------------------------------------------
	static async ls(args) {
		// Ignore any options
		const path = args.filter(arg => !arg.startsWith("-"))[0] || ".";
		const output = await CoreUtils.CLI.ls(path);

		// If path doesn't exist
		if(!output)
			return `${path}: No such file or directory`;

        // If the path isn't for a single file
		if(!output.mode)
            return output.join("\n");
	}
}
