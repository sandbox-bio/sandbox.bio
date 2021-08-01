// 	// -------------------------------------------------------------------------
// 	// ls
// 	// -------------------------------------------------------------------------
// 	static async ls(args, raw=false) {
// 		// Ignore flags (use try/catch b/c it's ok if user doesn't specify an arg with ls)
// 		let path;
// 		try {
// 			path = parseArgs(args)._[0];
// 		} catch (error) {}
// 		path = path || ".";

// 		// Get info about files in that path
// 		let stats = [];
// 		try {
// 			let stat = await CoreUtils.FS.stat(path);

// 			// If the path is a file, we already have the info we need
// 			if(!await CoreUtils.FS.isDir(stat.mode))
// 				stats = [{
// 					name: path.split("/").pop(),
// 					size: stat.size,
// 					date: stat.mtime.toLocaleString()
// 				}];

// 			// But if the path is a folder, get stat for each node in the folder
// 			else
// 			{
// 				const files = await CoreUtils.FS.readdir(path);
// 				for(let f of files) {
// 					if(f == "." || f == "..")
// 						continue;
// 					stat = await CoreUtils.FS.stat(`${path}/${f}`);
// 					const isDir = await CoreUtils.FS.isDir(stat.mode);
// 					let name = isDir ? f + "/" : f;
// 					name = (!raw && isDir) ? `\x1b[1;34m${name}\x1b[0m` : name;
// 					stats.push({
// 						name: name,
// 						size: stat.size,
// 						date: stat.mtime.toLocaleString()
// 					});
// 				}
// 			}

// 			if(raw)
// 				return stats;

// 			// Columnify output
// 			return columnify(stats.map(f => ({
// 				...f,
// 				size: prettyBytes(f.size),
// 				date: f.date.toLocaleString()
// 			})), {
// 				showHeaders: false,
// 				columnSplitter: "\t",
// 				columns: ["size", "date", "name"]
// 			});
// 		} catch (error) {
// 			console.error(error)
// 			throw `${path}: No such file or directory`;
// 		}
// 	}

// 	// Alias for ls :)
// 	static async ll() {
// 		return CoreUtils.ls(...arguments);
// 	}

// 	// -------------------------------------------------------------------------
// 	// Wrangling
// 	// -------------------------------------------------------------------------

// 	// grep: assume format = "grep <patternHere> <file>"
// 	static async grep(args)
// 	{
// 		args = parseArgs(args, { boolean: [ "E" ] });
// 		const pattern = args._[0].replaceAll('"', '').replaceAll("'", "");
// 		const path = args._[1];
// 		const contents = await CoreUtils.cat([path]);
// 		return contents.split("\n").filter(line => line.match(pattern)).join("\n");
// 	}

// 	// -------------------------------------------------------------------------
// 	// File system commands with external interaction
// 	// -------------------------------------------------------------------------

// 	// Generate a URL where you can download a file from the virtual file system
// 	static async download(args) {
// 		// Create blob from file contents
// 		const blob = new Blob([ await CoreUtils.FS.readFile(args[0]) ]);
// 		const url = URL.createObjectURL(blob);

// 		// Create a temp element to download it (changing window.location doesn't work for text files)
// 		let tmpLink = document.createElement("a");
// 		tmpLink.href = url;
// 		tmpLink.download = args[0];
// 		// Add it to the DOM and click on it
// 		document.body.appendChild(tmpLink);
// 		tmpLink.click();
// 		// Cleanup
// 		window.URL.revokeObjectURL(url);
// 		tmpLink.remove();

// 		return "Your download has started.";
// 	}
