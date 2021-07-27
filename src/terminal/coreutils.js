// =============================================================================
// Coreutils
// =============================================================================

import columnify from "columnify";
import prettyBytes from "pretty-bytes";

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
	static async ls(args, raw=false) {
		// Ignore flags
		const path = args.filter(arg => !arg.startsWith("-"))[0] || ".";

		// Get info about files in that path
		let stats = [];
		try {
			let stat = await CoreUtils.FS.stat(path);

			// If the path is a file, we already have the info we need
			if(!await CoreUtils.FS.isDir(stat.mode))
				stats = [{
					name: path.split("/").pop(),
					size: stat.size,
					date: stat.mtime.toLocaleString()
				}];

			// But if the path is a folder, get stat for each node in the folder
			else
			{
				const files = await CoreUtils.FS.readdir(path);
				for(let f of files) {
					if(f == "." || f == "..")
						continue;
					stat = await CoreUtils.FS.stat(`${path}/${f}`);
					const isDir = await CoreUtils.FS.isDir(stat.mode);
					let name = isDir ? f + "/" : f;
					name = (!raw && isDir) ? `\x1b[1;34m${name}\x1b[0m` : name;
					stats.push({
						name: name,
						size: stat.size,
						date: stat.mtime.toLocaleString()
					});
				}
			}

			if(raw)
				return stats;

			// Columnify output
			return columnify(stats.map(f => ({
				...f,
				size: prettyBytes(f.size),
				date: f.date.toLocaleString()
			})), {
				showHeaders: false,
				columnSplitter: "\t",
				columns: ["size", "date", "name"]
			});
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
	// Small utilities
	// -------------------------------------------------------------------------

	static async cat(args) {
		return await CoreUtils.CLI.cat(args[0]);
	}

	static async pwd() {
		return await CoreUtils.FS.cwd();
	}

	static async cd(args) {
		return await CoreUtils.FS.chdir(args[0]) || args[0];
	}

	static async echo(args) {
		return args.join(" ");
	}

	static async mv(args) {
		await CoreUtils.FS.rename(args[0], args[1]);
		return "";
	}

	static async mkdir(args) {
		// Don't use await because otherwise get "object could not be cloned" error
		args.map(d => CoreUtils.FS.mkdir(d));
		return "";
	}

	static async rmdir(args) {
		// Don't use await because otherwise get "object could not be cloned" error
		args.map(d => CoreUtils.FS.rmdir(d));
		return "";
	}
}
