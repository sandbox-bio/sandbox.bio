// =============================================================================
// Coreutils
// =============================================================================

// Simulating a few basic utils
export class CoreUtils
{
	static CLI;
	static get FS() {
		return CoreUtils.CLI.tools[1].module.FS;
	}

	// -------------------------------------------------------------------------
	// ls
	// -------------------------------------------------------------------------
	static async ls(args) {
		// Ignore flags
		const path = args.filter(arg => !arg.startsWith("-"))[0] || ".";

		// Get info about files in that path
		let stats = [];
		try {
			const info = await CoreUtils.FS.stat(path);

			// If the path is a file, we already have the info we need
			if(!await CoreUtils.FS.isDir(info.mode))
				stats = [{
					name: path.split("/").pop(),
					...info
				}];
			// But if the path is a folder, get stat for each node in the folder
			else
			{
				const files = await CoreUtils.FS.readdir(path);
				for(let f of files) {
					if(f == "." || f == "..")
						continue;
					stats.push({
						name: f,
						...(await CoreUtils.FS.stat(`${path}/${f}`))
					});
				}
			}

			return stats.map(f => {
				return `${f.name}\t${f.size}\t${f.mtime}`;
			}).join("\n");
		} catch (error) {
			console.error(error)
			throw `${path}: No such file or directory`;
		}
	}

	// Alias for ls :)
	static async ll() {
		return CoreUtils.ls(...arguments);
	}

	// -------------------------------------------------------------------------
	//
	// -------------------------------------------------------------------------
	static async pwd() {
		return await CoreUtils.FS.cwd();
	}

	static async chdir(args) {
		// Will return an error string if it doesn't work
		return await CoreUtils.FS.chdir(args[0]) || args[0];
	}
}
