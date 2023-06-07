const COREUTILS = ["ls", "head", "tail", "wc", "cat", "echo", "cut"];

export const TOOLS = [
	"samtools/1.10",
	{ loading: "eager", tool: "grep", version: "3.7", reinit: true  },
	...COREUTILS.map(program => ({ program, tool: "coreutils", loading: "eager", version: "8.32", reinit: true }))
]
