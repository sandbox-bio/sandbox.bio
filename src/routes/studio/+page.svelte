<script>
import * as marked from "marked";
import autosize from "svelte-autosize";
import autosize2 from "autosize";
import { onMount } from "svelte";
import { Button, Icon } from "sveltestrap";
import { LocalState, STATE_STUDIO } from "$src/utils";

// -----------------------------------------------------------------------------
// State
// -----------------------------------------------------------------------------

const updateState = (txt) => LocalState.set(STATE_STUDIO, txt);
$: updateState(tutorial);

// Initial values
let tutorial = "Markdown from a **tutorial** step goes here.";
let example = `<Alert>This is a **note** about this tutorial</Alert>

Here is a link: [Link Text](https://google.com)

Get a file listing:

<Execute command="ls -lah" />

Here's a long command:

<Execute command={"bedtools merge \\ -i exons.bed -o count,collapse \\ -d 90 -c 1,4"} />
`;

// Set state
onMount(async () => {
	const state = await LocalState.get(STATE_STUDIO);
	if (!state) return;
	tutorial = state;

	// Resize textarea manually
	setTimeout(() => {
		autosize2.update(document.querySelector("textarea"));
	}, 100);
});

// Download Markdown file
function download() {
	const file = new File([tutorial], "Step.md", { type: "text/plain" });
	const link = document.createElement("a");
	const url = URL.createObjectURL(file);

	link.href = url;
	link.download = file.name;
	document.body.appendChild(link);
	link.click();

	document.body.removeChild(link);
	window.URL.revokeObjectURL(url);
}

// -----------------------------------------------------------------------------
// Marked.js settings
// -----------------------------------------------------------------------------

// Open URLs in new tab
const renderer = {
	link(href, title, text) {
		return `<a href="${href}" target="_blank">${text}</a>`;
	}
};

// Parse <Alert> tags
const tagAlert = {
	name: "alert",
	level: "block",
	start(src) {
		return src.match(/^<Alert>/)?.index;
	},
	tokenizer(src) {
		const rule = /^<Alert>[\n]*(.*)[\n]*<\/Alert>/;
		const match = rule.exec(src);
		if (match) {
			const token = {
				type: "alert",
				raw: match[0],
				text: match[1].trim(),
				tokens: []
			};
			this.lexer.inline(token.text, token.tokens);
			return token;
		}
	},
	renderer(token) {
		return `
			<div class="alert alert-primary" role="alert">
				${this.parser.parseInline(token.tokens)}
			</div>
		`;
	}
};

// Parse <Execute> tags
// FIXME: inline doesn't work unless it's on the same line (removing ^ doesn't quite work as intended)
const tagExecute = {
	name: "execute",
	level: "block",
	start(src) {
		return src.match(/^<Execute/)?.index;
	},
	tokenizer(src) {
		// return false;
		const rule = /^<Execute command="(.+?)" (inline )?\/>/;
		const match = rule.exec(src);
		if (match) {
			const token = {
				type: "execute",
				raw: match[0],
				text: match[1].trim(),
				tokens: [],
				isInline: match[2]?.trim()
			};
			this.lexer.inline(token.text, token.tokens);
			return token;
		}
	},
	renderer(token) {
		// Here we don't parseInline b/c everything inside that component should be raw text
		const command = token.text.replace(/ \\ /g, " \\ <br />&nbsp;&nbsp;&nbsp;");
		if (token.isInline) {
			return `
				<kbd>${command}</kbd>
			`;
		} else {
			return `
				<div class="card bg-dark text-white font-monospace mb-3">
					<div class="card-body">
						${command}
					</div>
				</div>
			`;
		}
	}
};

marked.use({ renderer, extensions: [tagAlert, tagExecute] });
</script>

<svelte:head>
	<title>Studio - sandbox.bio</title>
</svelte:head>

<!-- Tutorial -->
<div class="row">
	<div class="col-6">
		<h4>
			Markdown
			<Button color="primary" size="sm" on:click={download}>
				<Icon name="download" /> Download
			</Button>
		</h4>
		<textarea use:autosize class="form-control border border-primary" bind:value={tutorial} />
	</div>
	<div class="col-6">
		<h4>Preview</h4>
		<div class="border border-primary p-2">
			{@html marked.parse(tutorial)}
		</div>
	</div>
</div>

<hr />

<div class="row mt-4 opacity-75">
	<div class="col-6">
		<h4>Examples:</h4>
		<textarea use:autosize class="form-control border border-secondary" bind:value={example} />
	</div>
	<div class="col-6">
		<h4>Preview</h4>
		<div class="border border-secondary p-2">
			{@html marked.parse(example)}
		</div>
	</div>
</div>
